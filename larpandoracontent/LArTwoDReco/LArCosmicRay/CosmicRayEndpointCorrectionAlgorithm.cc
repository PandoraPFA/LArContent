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
    m_maxSmoothPoints(2),
    m_thresholdMaxAngleDeviation(25.f)
{
}

StatusCode CosmicRayEndpointCorrectionAlgorithm::Run()
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
    
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    while (!clusterVector.empty())
    {
        const Cluster *pCluster(clusterVector.front());
        
        TwoDSlidingFitResultMap::const_iterator microFitResult(microSlidingFitResultMap.find(pCluster));
        TwoDSlidingFitResultMap::const_iterator macroFitResult(macroSlidingFitResultMap.find(pCluster));

        if (microFitResult == microSlidingFitResultMap.end() || macroFitResult == macroSlidingFitResultMap.end())
            continue;        

        this->InvestigateClusterEndpoint(pCluster, true, microFitResult->second, macroFitResult->second);
        this->InvestigateClusterEndpoint(pCluster, false, microFitResult->second, macroFitResult->second);

        // THEN DELETE CLUSTER OUT OF MAP AND SLIDING FIT RESULT
        UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
    }
        
    return STATUS_CODE_SUCCESS;
}    
    
//------------------------------------------------------------------------------------------------------------------------------------------    

    
void CosmicRayEndpointCorrectionAlgorithm::InvestigateClusterEndpoint(const Cluster *const pCluster, const bool isUpstream, const TwoDSlidingFitResult &microFitResult, const TwoDSlidingFitResult &macroFitResult)
{   
    CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
    if (!GetClusterMergingCoordinates(microFitResult, macroFitResult, macroFitResult, isUpstream, clusterMergePoint, clusterMergeDirection))
    {
        std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return;
    }
    
    const CartesianVector &endpointPosition(isUpstream ? microFitResult.GetGlobalMaxLayerPosition() : microFitResult.GetGlobalMinLayerPosition());
    const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX()), predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
    const CartesianVector extrapolatedPoint(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));

    const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());
    std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;

    if (endpointSeparation < std::numeric_limits<float>::epsilon())
    {
        std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return;
    }

    ClusterAssociation clusterAssociation(isUpstream ? ClusterAssociation(pCluster, pCluster, clusterMergePoint, clusterMergeDirection, extrapolatedPoint, clusterMergeDirection * (-1.f)) :
                                                       ClusterAssociation(pCluster, pCluster, extrapolatedPoint, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f)));

    // DOES ENDPOINT NEED REMOVING?
    if (!this->IsCosmicRay(clusterAssociation, endpointPosition, isUpstream))
    {
        std::cout << "DOESN'T SATISY CRITERIA" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());              
        return;
    }

    // REMOVE HITS - CREATE AN ACTUAL CLUSTER ASSOCIATION HERE
    std::cout << "REMOVING HITS..." << std::endl;

    ClusterToCaloHitListMap clusterToCaloHitListMap;    
    GetExtrapolatedCaloHits(clusterAssociation, pClusterList, clusterToCaloHitListMap);

    // VISUALIZE COLLECTED CALO HITS
    for (auto entry : clusterToCaloHitListMap)
    {
        const CaloHitList &caloHitList(entry.second);

        for (auto &hit :caloHitList}
        {
            const CartesianVector &hitPosition(hit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POSITION", GREEN, 2);
        }
    }
            
    PandoraMonitoringApi::ViewEvent(this->GetPandora());


}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void TrackInEMShowerAlgorithm::RefineTrackEndpoint(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList *const pClusterList, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterVector &clusterVector) const
{
    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (const auto &entry : clusterToCaloHitListMap)
    {
        if (entry.first != clusterAssociation.GetUpstreamCluster())
            showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ClusterList remnantClusterList;

    CartesianVector &clusterMergePoint(isUpstream ? clusterAssociation.GetUpstreamMergePoint() : clusterAssociation.GetDownstreamMergePoint());
    const Cluster *const pMainTrackCluster(RemoveOffAxisHitsFromTrack(pCluster, clusterMergePoint, isUpstream, clusterToCaloHitListMap, remnantClusterList,
        microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));

    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        
        this->UpdateForClusterDeletion(pShowerCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        this->AddHitsToMainTrack(pShowerCluster, pMainTrackCluster, caloHitsToMerge, clusterAssociation, remnantClusterList);
    }

    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList);
    this->UpdateAfterMainTrackCreation(pMainTrackCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------ 


//------------------------------------------------------------------------------------------------------------------------------------------    

void CosmicRayEndpointCorrectionAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------


bool CosmicRayEndpointCorrectionAlgorithm::IsCosmicRay(const ClusterAssociation &clusterAssociation, const CartesianVector &clusterEndpoint, const bool isUpstream) const
{
    const CartesianVector &clusterMergePoint(isUpstream ? clusterAssociation.GetUpstreamMergePoint() : clusterAssociation.GetDownstreamMergePoint());
    const CartesianVector &clusterMergeDirection(isUpstream ? clusterAssociation.GetUpstreamMergeDirection() : clusterAssociation.GetDownstreamMergeDirection());
    const CartesianVector &clusterExtrapolatedPoint(isUpstream ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());

    //////////////////////////
    ClusterList theCluster({clusterAssociation.GetUpstreamCluster()});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterEndpoint, "ENDPOINT", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "CLUSTER MERGE POINT", RED, 2);
    //////////////////////////
       
   
    // INVESTIGATE ENDPOINT Z SEPARATION
    const float deltaZ(std::fabs(clusterEndpoint.GetZ() - clusterExtrapolatedPoint.GetZ()));
    const float extrapolatedSeparation((clusterExtrapolatedPoint - clusterMergePoint).GetMagnitude());
    const float scaledDeltaZ(deltaZ / extrapolatedSeparation); 
    
    /////////////////
    const CartesianVector start(clusterMergePoint + (clusterMergeDirection*20));
    const CartesianVector end(clusterMergePoint - (clusterMergeDirection*20));
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterExtrapolatedPoint, "JAM", GREEN, 2);
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

    const bool isClusterEndpointInBoundary((clusterEndpoint.GetX() < allowanceLowXEdge) || (clusterEndpoint.GetX() > allowanceHighXEdge));
    const bool isClusterMergePointInBoundary((clusterMergePoint.GetX() < allowanceLowXEdge) || (clusterMergePoint.GetX() > allowanceHighXEdge));

    std::cout << "Distance from boundary: " << std::min(std::fabs(clusterMergePoint.GetX() - tpcHighXEdge), std::fabs(clusterMergePoint.GetX() - tpcLowXEdge)) << std::endl;

    if (!(isClusterEndpointInBoundary || isClusterMergePointInBoundary))
    {
        std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

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

    bool reachedFirstCurve(false), isClockwise(false), metCriteria(false);
    float previousOpeningAngle(std::numeric_limits<float>::max());
    unsigned int tooSmooth(0);
    
    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        float microOpeningAngle(microDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);

        const bool isAbove(microDirection.GetZ() > (clusterAssociation.GetConnectingLineDirection().GetZ() * microDirection.GetX() / clusterAssociation.GetConnectingLineDirection().GetX()));
        if (!isAbove)
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
        if (std::fabs(microOpeningAngle > m_thresholdAngleDeviation) && !reachedFirstCurve)
        {
            reachedFirstCurve = true;
            isClockwise = (microOpeningAngle > 0.f);
            continue;
        }

        if (reachedFirstCurve)
        {
            if ((isClockwise && (microOpeningAngle < previousOpeningAngle)) || (!isClockwise && (microOpeningAngle > previousOpeningAngle)) ||
                (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
            {
                ++tooSmooth;

                if (tooSmooth > m_maxSmoothPoints)
                {
                    std::cout << "TOO SMOOTH" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    return false;
                }
            }
            else
            {
                tooSmooth = 0;
            
                if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
                {
                    std::cout << "MEET ANGLE CRITERIA" << std::endl;
                    metCriteria = true;
                    break;
                }
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
    
    return true;

}

//------------------------------------------------------------------------------------------------------------------------------------------    


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




StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromTPC", m_maxDistanceFromTPC));

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
        "MaxSmoothPoints", m_maxSmoothPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdMaxAngleDeviation", m_thresholdMaxAngleDeviation));           

        
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content


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
