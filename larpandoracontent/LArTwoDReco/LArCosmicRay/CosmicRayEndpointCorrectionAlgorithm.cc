/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.cc 
 *
 *  @brief  Implementation of the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{
    
CosmicRayEndpointCorrectionAlgorithm::CosmicRayEndpointCorrectionAlgorithm() :
    m_minCaloHits(50),
    m_maxDistanceFromTPC(15.f),
    m_curveThreshold(0.3f),
    m_minScaledZOffset(0.25f),
    m_thresholdAngleDeviation(10.f),
    m_thresholdAngleDeviationBetweenLayers(1.f),
    m_maxAnomalousPoints(2),
    m_thresholdMaxAngleDeviation(25.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------        

StatusCode CosmicRayEndpointCorrectionAlgorithm::Run()
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    ClusterVector clusterVector;
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    SlidingFitResultMapPair slidingFitResultMapVector({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    
    this->InitialiseContainers(pClusterList, clusterVector, slidingFitResultMapVector);

    while (!clusterVector.empty())
    {
        const Cluster *pCluster(clusterVector.front());

        //COULD CHANGE TO AN AT
        TwoDSlidingFitResultMap::const_iterator microFitResult(slidingFitResultMapVector.first->find(pCluster));
        TwoDSlidingFitResultMap::const_iterator macroFitResult(slidingFitResultMapVector.second->find(pCluster));

        if (microFitResult == microSlidingFitResultMap.end() || macroFitResult == macroSlidingFitResultMap.end())
        {
            std::cout << "THIS SHOULD NEVER EVER EVER HAPPEN" << std::endl;
            continue;
        }

        ClusterToCaloHitListMap clusterToCaloHitListMap;
        
        ClusterAssociation upstreamClusterAssociation;
        const bool modifyUpstreamEnd(this->InvestigateClusterEnd(pCluster, false, microFitResult->second, macroFitResult->second, upstreamClusterAssociation));

        if (modifyUpstreamEnd)
            GetExtrapolatedCaloHits(upstreamClusterAssociation, pClusterList, clusterToCaloHitListMap);

        ClusterAssociation downstreamClusterAssociation;
        const bool modifyDownstreamEnd(this->InvestigateClusterEnd(pCluster, true, microFitResult->second, macroFitResult->second, upstreamClusterAssociation));

        if (modifyDownstreamEnd)
            GetExtrapolatedCaloHits(downstreamClusterAssociation, pClusterList, clusterToCaloHitListMap);

        
        if (modifyUpstreamEnd || modifyDownstreamEnd)
            this->RefineTrackEndpoint(pCluster, downstreamClusterAssociation.GetDownstreamMergePoint(), upstreamClusterAssociation.GetUpstreamMergePoint(),
                clusterToCaloHitListMap, pClusterList, clusterVector, slidingFitResultMapVector);
    }
        
    return STATUS_CODE_SUCCESS;
}    
    
//------------------------------------------------------------------------------------------------------------------------------------------    
    
bool CosmicRayEndpointCorrectionAlgorithm::InvestigateClusterEnd(const Cluster *const pCluster, const bool isUpstream, const TwoDSlidingFitResult &microFitResult, const TwoDSlidingFitResult &macroFitResult, ClusterAssociation &clusterAssociation)
{   
    CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
    if (!GetClusterMergingCoordinates(microFitResult, macroFitResult, macroFitResult, isUpstream, clusterMergePoint, clusterMergeDirection))
    {
        std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }
    
    const CartesianVector &endpointPosition(isUpstream ? microFitResult.GetGlobalMaxLayerPosition() : microFitResult.GetGlobalMinLayerPosition());
    const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX()), predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
    const CartesianVector extrapolatedEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));

    const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());
    std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;

    if (endpointSeparation < std::numeric_limits<float>::epsilon())
    {
        std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    // INVESTIGATE ENDPOINT Z SEPARATION
    const float deltaZ(std::fabs(endpointPosition.GetZ() - extrapolatedEndpointPosition.GetZ()));
    const float extrapolatedSeparation((extrapolatedEndpointPosition - clusterMergePoint).GetMagnitude());
    const float scaledDeltaZ(deltaZ / extrapolatedSeparation); 
    
    /////////////////
    const CartesianVector start(clusterMergePoint + (clusterMergeDirection*20));
    const CartesianVector end(clusterMergePoint - (clusterMergeDirection*20));
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedEndpointPosition, "JAM", GREEN, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "JAM LINE", GREEN, 2, 2);
    std::cout << "deltaZ: " << deltaZ << std::endl;
    std::cout << "scaled deltaZ: " << scaledDeltaZ << std::endl;    
    /////////////////

    if (scaledDeltaZ < m_minScaledZOffset)
    {
        std::cout << "SCALED DELTA Z IS TOO HIGH" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    // INVESTIGATE PROXIMITY TO TPC BOUNDARY - IGNORE CASES WHERE IT IS NEAR THE EXTREME EDGE (ISOBEL TODO)
    const LArTPC &pLArTPC(this->GetPandora().GetGeometry()->GetLArTPC());
    const float tpcHighXEdge(pLArTPC.GetCenterX() + (pLArTPC.GetWidthX() / 2.f)), tpcLowXEdge(pLArTPC.GetCenterX() - (pLArTPC.GetWidthX() / 2.f));
    const float allowanceHighXEdge(tpcHighXEdge - m_maxDistanceFromTPC), allowanceLowXEdge(tpcLowXEdge + m_maxDistanceFromTPC);

    const bool isClusterEndpointInBoundary((endpointPosition.GetX() < allowanceLowXEdge) || (endpointPosition.GetX() > allowanceHighXEdge));
    const bool isClusterMergePointInBoundary((clusterMergePoint.GetX() < allowanceLowXEdge) || (clusterMergePoint.GetX() > allowanceHighXEdge));

    std::cout << "Distance from boundary: " << std::min(std::fabs(clusterMergePoint.GetX() - tpcHighXEdge), std::fabs(clusterMergePoint.GetX() - tpcLowXEdge)) << std::endl;

    if (!(isClusterEndpointInBoundary || isClusterMergePointInBoundary))
    {
        std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    clusterAssociation = isUpstream ? ClusterAssociation(pCluster, pCluster, clusterMergePoint, clusterMergeDirection, extrapolatedEndpointPosition, clusterMergeDirection * (-1.f)) :
                                                       ClusterAssociation(pCluster, pCluster, extrapolatedEndpointPosition, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f));
    
    // DOES ENDPOINT NEED REMOVING?
    if(!this->IsDeltaRay(clusterAssociation, isUpstream))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayEndpointCorrectionAlgorithm::IsDeltaRay(const ClusterAssociation &clusterAssociation, const bool isUpstream) const
{
    const CartesianVector &clusterMergePoint(isUpstream ? clusterAssociation.GetUpstreamMergePoint() : clusterAssociation.GetDownstreamMergePoint());
    const CartesianVector &clusterMergeDirection(isUpstream ? clusterAssociation.GetUpstreamMergeDirection() : clusterAssociation.GetDownstreamMergeDirection());

    //////////////////////////
    ClusterList theCluster({clusterAssociation.GetUpstreamCluster()});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "CLUSTER MERGE POINT", RED, 2);
    //////////////////////////
       
    // MAKE MORE DETAILED FIT OF THE AMBIGUOUS SECTION
    CartesianPointVector hitSubset;
    const OrderedCaloHitList &orderedCaloHitList(clusterAssociation.GetUpstreamCluster()->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if ((isUpstream && (hitPosition.GetZ() > clusterMergePoint.GetZ())) || (!isUpstream && (hitPosition.GetZ() < clusterMergePoint.GetZ())))
                hitSubset.push_back(hitPosition);          
        }
    }

    TwoDSlidingFitResultList subsetFitVector;
    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const int m_slidingFitWindow(5);
        const TwoDSlidingFitResult subsetFit(&hitSubset, m_slidingFitWindow, slidingFitPitch);

        subsetFitVector.push_back(subsetFit);
    }
    catch (const StatusCodeException &)
    {
        std::cout << "CANNOT MAKE A FIT" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    // INVESTIGATE THE DIRECTION CHANGE
    const TwoDSlidingFitResult subsetFit(subsetFitVector.front());
    const LayerFitResultMap &clusterMicroLayerFitResultMap(subsetFit.GetLayerFitResultMap());
    const int startLayer(isUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    const int endLayer(isUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    const int loopTerminationLayer(isUpstream ? endLayer + 1 : endLayer - 1);
    const int step(isUpstream ? 1 : -1);

    unsigned int anomalousLayerCount(0);    
    bool reachedFirstCurve(false), isCurveClockwise(false);
    float previousOpeningAngle(std::numeric_limits<float>::max());
    
    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        float microOpeningAngle(microDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);

        if(microDirection.GetZ() < (clusterAssociation.GetConnectingLineDirection().GetZ() * microDirection.GetX() / clusterAssociation.GetConnectingLineDirection().GetX()))
            microOpeningAngle *= (-1.f);

        const float layerAngleDeviation(previousOpeningAngle > 180.f ? microOpeningAngle : std::fabs(microOpeningAngle - previousOpeningAngle));
        
        /////////////////////
        CartesianVector microPosition(0.f, 0.f, 0.f);
        subsetFit.GetGlobalPosition(microIter->second.GetL(), microIter->second.GetFitT(), microPosition);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &microPosition, "MICRO POSITION", BLACK, 2);
        std::cout << "ANGLE DEVIATION: " << microOpeningAngle << std::endl;
        std::cout << "<-------------------- CHANGE FROM LAST LAYER: " << layerAngleDeviation << std::endl;
        /////////////////////
        
        // ISOBEL - DO YOU NEED FINAL THRESHOLD
        if ((!reachedFirstCurve) && (std::fabs(microOpeningAngle) > m_thresholdAngleDeviation))
        {
            reachedFirstCurve = true;
            isCurveClockwise = (microOpeningAngle > 0.f);
            continue;
        }

        if (reachedFirstCurve)
        {
            if ((isCurveClockwise && (microOpeningAngle < previousOpeningAngle)) || (!isCurveClockwise && (microOpeningAngle > previousOpeningAngle)) ||
                (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
            {
                ++anomalousLayerCount;

                if (anomalousLayerCount > m_maxAnomalousPoints)
                {
                    std::cout << "TOO SMOOTH" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    return false;
                }
            }
            else
            {
                anomalousLayerCount = 0;
            
                if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
                {
                    std::cout << "MEET ANGLE CRITERIA" << std::endl;
                    return true;
                }
            }
        }

        previousOpeningAngle = microOpeningAngle;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void CosmicRayEndpointCorrectionAlgorithm::RefineTrackEndpoint(const Cluster *const pCluster, const CartesianVector &clusterUpstreamMergePoint, const CartesianVector &clusterDownstreamMergePoint, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (const auto &entry : clusterToCaloHitListMap)
    {
        if (entry.first != pCluster)
            showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ClusterList remnantClusterList;

    const Cluster *const pIntermediateCluster(RemoveOffAxisHitsFromTrack(pCluster, clusterUpstreamMergePoint, false, clusterToCaloHitListMap, remnantClusterList,
        *slidingFitResultMapPair.first, *slidingFitResultMapPair.second));
                                              
    const Cluster *const pFinalCluster(RemoveOffAxisHitsFromTrack(pCluster, clusterDownstreamMergePoint, true, clusterToCaloHitListMap, remnantClusterList,
        *slidingFitResultMapPair.first, *slidingFitResultMapPair.second));

    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        this->AddHitsToMainTrack(pShowerCluster, pFinalCluster, caloHitsToMerge, clusterAssociation, remnantClusterList);
    }

    ClusterList createdClusters;
    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList, createdClusters);

    showerClustersToFragment.push_back(pCluster);
    this->UpdateContainers(showerClustersToFragment, createdClusters, clusterVector, slidingFitResultMapPair);
}
   
//------------------------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromTPC", m_maxDistanceFromTPC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurveThreshold", m_curveThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinScaledZOffset", m_minScaledZOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdAngleDeviation", m_thresholdAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdAngleDeviationBetweenLayers", m_thresholdAngleDeviationBetweenLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAnomalousPoints", m_maxAnomalousPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdMaxAngleDeviation", m_thresholdMaxAngleDeviation));           

    return CosmicRayTrackRefinementBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content



/*
    // INVESTIGATE THE STRAIGHTNESS OF THE CLUSTER SUBSET
    const CartesianVector &minPosition(subsetFitVector.front().GetGlobalMinLayerPosition()), &maxPosition(subsetFitVector.front().GetGlobalMaxLayerPosition());
    CartesianVector straightLinePrediction(maxPosition - minPosition);
    straightLinePrediction = straightLinePrediction.GetUnitVector();
    
    float distanceFromLine(0.f);
    unsigned int distanceCount(0);
    for (const CartesianVector &hitPosition : hitSubset)
    {
        if (hitPosition.GetZ() < minPosition.GetZ() || hitPosition.GetZ() > maxPosition.GetZ())
            continue;

        distanceFromLine += (straightLinePrediction.GetCrossProduct(hitPosition - minPosition).GetMagnitude());

        ++distanceCount;
    }

    if (distanceCount == 0)
    {
        std::cout << "DISTANCE COUNT TOO SMALL" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    const float averageDistanceFromLine(distanceFromLine / distanceCount);
    std::cout << "averageDistanceFromLine: " << averageDistanceFromLine << std::endl;    
    
    if (averageDistanceFromLine < m_curveThreshold)
    {
        std::cout << "TOO STRAIGHT" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }
*/




///OLD METHOD
/*

        bool isAbove(microDirection.GetZ() > (predictedGradient*microDirection.GetX()));
        if (!isAbove)
            microOpeningAngle *= (-1.f);

        std::cout << "ANGLE DEVIATION: " << microOpeningAngle << std::endl;

        float layerAngleDeviation(0.f);
        if (reachedFirstCurve)
        {
            layerAngleDeviation = (std::fabs(microOpeningAngle - previousOpeningAngle));
            std::cout << "<-------------------- CHANGE FROM LAST LAYER: " << layerAngleDeviation << std::endl;
        }
        
        if (((microOpeningAngle > m_thresholdAngleDeviation) || (microOpeningAngle < -m_thresholdAngleDeviation)) && !reachedFirstCurve)
        {
            reachedFirstCurve = true;
            isClockwise = (microOpeningAngle < 0.f);
            previousOpeningAngle = microOpeningAngle;
            continue;
        }

        if (reachedFirstCurve)
        {
            if (((isClockwise && (microOpeningAngle > previousOpeningAngle)) || (!isClockwise && (microOpeningAngle < previousOpeningAngle))) && (layerAngleDeviation > 1 ))
            {
                ++numberOfWobbles;

                if (numberOfWobbles > m_maxNumberOfWobbles)
                {
                    std::cout << "TOO MANY WOBBLES - RESET" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    return false;
                }
            }

            if ((isClockwise && (microOpeningAngle < -25.f)) || (!isClockwise && (microOpeningAngle > 25.f)))
            {
                std::cout << "MEET ANGLE CRITERIA" << std::endl;
                metCriteria = true;
                break;
            }
        }

        previousOpeningAngle = microOpeningAngle;
    }


    if(!metCriteria)
    {
        std::cout << "DIDN'T FIND ANGLE CRITERIA" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }



 */
