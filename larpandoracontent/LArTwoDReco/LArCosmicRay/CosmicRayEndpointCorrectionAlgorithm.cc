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

/*
        const Cluster *const pClusterToRemove(clusterAssociationCaloHitOwnershipMap.begin()->first.GetUpstreamCluster()->GetNCaloHits() >
                clusterAssociationCaloHitOwnershipMap.begin()->first.GetDownstreamCluster()->GetNCaloHits() ?
                clusterAssociationCaloHitOwnershipMap.begin()->first.GetDownstreamCluster() :
                clusterAssociationCaloHitOwnershipMap.begin()->first.GetUpstreamCluster());
 */

    
StatusCode CosmicRayEndpointCorrectionAlgorithm::Run()
{
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    ClusterVector clusterVector;
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    SlidingFitResultMapPair slidingFitResultMapPair({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    
    this->InitialiseContainers(pClusterList, clusterVector, slidingFitResultMapPair);

    //CONSIDER CLUSTERS, ONCE CONSIDERED REMOVE FROM THE CLUSTER VECTOR 
    
    bool mergeMade(false);
    unsigned int loopIterations(0);
    do
    {
        ++loopIterations;

        ClusterAssociationVector clusterAssociationVector;
        this->FindBestClusterAssociation(clusterVector, slidingFitResultMapPair, clusterAssociationVector);
        
        if (clusterAssociationVector.empty())
            break;

        for (ClusterAssociation &clusterAssociation : clusterAssociationVector)
        {
            ClusterToCaloHitListMap clusterToCaloHitListMap;
            this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, clusterToCaloHitListMap);

            this->IsTrackContinuous(clusterAssociationCaloHitOwnershipMap);
        /* 
        if (clusterAssociationCaloHitOwnershipMap.empty())
        {
            RemoveClusterFromContainers(pClusterToRemove, clusterVector, slidingFitResultMapPair);
            continue;
        }
        */
        //NEED TO CHECK IF CLUSTERS ARE CONTAINED IN BOTH LISTS.
        
        /*
        this->CreateMainTrack(clusterAssociation, clusterToCaloHitListMap, pClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector);

        mergeMade = true;
        */
        }
    }
    while(mergeMade && (loopIterations < 10));

    return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------    
    
void CosmicRayEndpointCorrectionAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    ClusterAssociationVector &clusterAssociationVector)
{
    //ATTN: Method assumes clusters ordered by their hits 
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));
        
        for (unsigned int isEndUpstream = 0; isEndUpstream < 2; ++isEndUpstream)
        {
            const bool isClusterUpstream(!isEndUpstream);
            
            CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
            if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isClusterUpstream, clusterMergePoint, clusterMergeDirection))
            {
                std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }
    
            const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
            const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX()), predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
            const CartesianVector extrapolatedEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));

            const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());
            std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;

            if (endpointSeparation < std::numeric_limits<float>::epsilon())
            {
                std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // INVESTIGATE ENDPOINT Z SEPARATION
            const float deltaZ(std::fabs(endpointPosition.GetZ() - extrapolatedEndpointPosition.GetZ()));
            const float extrapolatedSeparation((extrapolatedEndpointPosition - clusterMergePoint).GetMagnitude());
            const float scaledDeltaZ(deltaZ / extrapolatedSeparation);;
    
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
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
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
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // DOES ENDPOINT NEED REMOVING?
            
            if(this->IsDeltaRay(pCluster, clusterMergePoint, clusterMergeDirection, isEndUpstream))
            {
                ClusterEndpointAssociation clusterEndpointAssociation(isEndUpstream ?
                    ClusterEndpointAssociation(extrapolatedEndpointPosition, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
                    ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedEndpointPosition, clusterMergeDirection * (-1.f), pCluster, false));
                
                clusterAssociationVector.push_back(clusterEndpointAssociation);
            }
        }

        if (!clusterAssociationVector.empty())
            break;
                
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayEndpointCorrectionAlgorithm::IsDeltaRay(const Cluster *const pCluster, const CartesianVector &clusterMergePoint,
    const CartesianVector &clusterMergeDirection, const bool isEndUpstream) const
{
    // MAKE MORE DETAILED FIT OF THE AMBIGUOUS SECTION
    CartesianPointVector hitSubset;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if ((isEndUpstream && (hitPosition.GetZ() < clusterMergePoint.GetZ())) || (!isEndUpstream && (hitPosition.GetZ() > clusterMergePoint.GetZ())))
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
    const int startLayer(isEndUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    const int endLayer(isEndUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    const int loopTerminationLayer(isEndUpstream ? endLayer - 1 : endLayer + 1);
    const int step(isEndUpstream ? -1 : 1);

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

        // IT IS USING THE CORRECT DEFINITION - AT THIS POINT ALL DIRECTIONS HAVE POSITIVE Z?
        if(microDirection.GetZ() < (clusterMergeDirection.GetZ() * microDirection.GetX() / clusterMergeDirection.GetX()))
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
/*
void CosmicRayEndpointCorrectionAlgorithm::CreateMainTrack(const ClusterAssociationCaloHitOwnershipMap &clusterAssociationCaloHitOwnershipMap, const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    ClusterAssociationVector clusterAssociationVector;
    for (auto &entry : clusterAssociationCaloHitOwnershipMap)
        clusterAssociationVector.push_back(&entry.first);

    for (const ClusterAssociation *const clusterAssociation : clusterAssociationVector)
    {
        const ClusterToCaloHitListMap &clusterToCaloHitListMap(clusterAssociationCaloHitOwnershipMap.at(*clusterAssociation));
        
        // Determine the shower clusters which contain hits that belong to the main track
        ClusterVector showerClustersToFragment;
        for (auto &entry : clusterToCaloHitListMap)
        {
            if (entry.first != clusterAssociation->GetUpstreamCluster() || entry.first != clusterAssociation->GetDownstreamCluster())
                showerClustersToFragment.push_back(entry.first);
        }

        std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

        ClusterList remnantClusterList;

        const Cluster *pMainTrackCluster(clusterAssociation->GetUpstreamCluster());
        const TwoDSlidingFitResult &microFitResult(slidingFitResultMapPair.first->at(pMainTrackCluster));
        
        bool isUpstreamEnd(std::min(clusterAssociation->GetUpstreamMergePoint().GetDistanceSquared(microFitResult.GetGlobalMinLayerPosition()),
                                    clusterAssociation->GetDownstreamMergePoint().GetDistanceSquared(microFitResult.GetGlobalMinLayerPosition())) <
                           std::min(clusterAssociation->GetUpstreamMergePoint().GetDistanceSquared(microFitResult.GetGlobalMaxLayerPosition()),
                                    clusterAssociation->GetDownstreamMergePoint().GetDistanceSquared(microFitResult.GetGlobalMaxLayerPosition())));

        if (isUpstreamEnd)
        {
            pMainTrackCluster = RemoveOffAxisHitsFromTrack(pMainTrackCluster, clusterAssociation->GetDownstreamMergePoint(), isUpstreamEnd,
                clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second);
        }
        else
        {
            pMainTrackCluster = RemoveOffAxisHitsFromTrack(pMainTrackCluster, clusterAssociation->GetUpstreamMergePoint(), isUpstreamEnd,
                clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second);
        }

        for (const Cluster *const pShowerCluster : showerClustersToFragment)
        {
            const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
            this->AddHitsToMainTrack(pShowerCluster, pMainTrackCluster, caloHitsToMerge, *clusterAssociation, remnantClusterList);
        }
    
        ClusterList createdClusters;
        this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList, createdClusters);

        showerClustersToFragment.push_back(clusterAssociation->GetUpstreamCluster());
        this->UpdateContainers(showerClustersToFragment, createdClusters, clusterVector, slidingFitResultMapPair);

        //IF CV MORE THAN 2 NEED TO UPDATE AFTER MODIFICATION
    }
}
   
//------------------------------------------------------------------------------------------------------------------------------------------
void CosmicRayEndpointCorrectionAlgorithm::UpdateAfterMainModification(const Cluster *const pDeletedCluster, const Cluster *const pNewCluster, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterAssociationCaloHitOwnershipMap &clusterAssociationCaloHitOwnershipMap) const
{
    // REPLACE IN MICRO FIT MAP
    const TwoDSlidingFitResultMap::const_iterator microFitToDelete(slidingFitResultMapPair.first->find(pDeletedCluster));
    
    if (microFitToDelete == slidingFitResultMapPair.first->end())
    {
        std::cout << "ISOBEL THIS IS BAD" << std::endl;
        throw STATUS_CODE_FAILURE;
    }
    
    slidingFitResultMapPair.first->insert(std::make_pair(pNewCluster, microFitToDelete->second));
    slidingFitResultMapPair.first->erase(microFitToDelete);

    // REPLACE IN MACRO FIT MAP
    const TwoDSlidingFitResultMap::const_iterator macroFitToDelete(slidingFitResultMapPair.second->find(pDeletedCluster));
    
    if (macroFitToDelete == slidingFitResultMapPair.second->end())
    {
        std::cout << "ISOBEL THIS IS BAD" << std::endl;
        throw STATUS_CODE_FAILURE;
    }
    
    slidingFitResultMapPair.second->insert(std::make_pair(pNewCluster, macroFitToDelete->second));
    slidingFitResultMapPair.second->erase(microFitToDelete);

    ClusterAssociationVector clusterAssociationVector;
    for (auto &entry : clusterAssociationCaloHitOwnershipMap)
        clusterAssociationVector.push_back(&entry.first);
    
    for (const ClusterAssociation *clusterAssociation : clusterAssociationVector)
    {
        if (pDeletedCluster == clusterAssociation->GetUpstreamCluster())
            clusterAssociation->SetUpstreamCluster(pNewCluster);

        if (pDeletedCluster == clusterAssociation->GetDownstreamCluster())
            clusterAssociation->SetDownstreamCluster(pNewCluster);
        
        // then alter this
        //const ClusterToCaloHitListMap &clusterToCaloHitListMap(clusterAssociationCaloHitOwnershipMap.at(*clusterAssociation));
    }

}
*/


//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayEndpointCorrectionAlgorithm::ClusterEndpointAssociation::ClusterEndpointAssociation(const CartesianVector &upstreamMergePoint, const CartesianVector &upstreamMergeDirection,
    const CartesianVector &downstreamMergePoint, const CartesianVector &downstreamMergeDirection, const Cluster *pMainTrackCluster) :
        ClusterAssociation(upstreamMergePoint, upstreamMergeDirection, downstreamMergePoint, downstreamMergeDirection),
        m_pMainTrackCluster(pMainTrackCluster)
{
}

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
