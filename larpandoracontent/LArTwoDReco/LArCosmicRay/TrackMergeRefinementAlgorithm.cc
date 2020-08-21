/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackMergeRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the track refinement class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackMergeRefinementAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"

using namespace pandora;

namespace lar_content
{

TrackMergeRefinementAlgorithm::TrackMergeRefinementAlgorithm() :
    m_minClusterLengthSum(75.f),
    m_minSeparationDistance(0.f),
    m_minDirectionDeviationCosAngle(0.99f),
    m_maxPredictedMergePointOffset(5.f),
    m_distanceFromLine(0.35f),
    m_boundaryTolerance(2.f)  
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    

    
StatusCode TrackMergeRefinementAlgorithm::Run()
{

    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);    

    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    ClusterVector clusterVector;
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    SlidingFitResultMapPair slidingFitResultMapPair({&microSlidingFitResultMap, &macroSlidingFitResultMap});

    this->InitialiseContainers(pClusterList, LArClusterHelper::SortByNHits, clusterVector, slidingFitResultMapPair);

    unsigned int loopIterations(0);
    ClusterList createdMainTrackClusters;    
    while(loopIterations < 10) 
    {
        ++loopIterations;

        std::cout << "\033[31m" << "TrackMerge: Finding best cluster association..." << "\033[0m"  << std::endl;
        ClusterPairAssociation clusterAssociation;
        if(!this->FindBestClusterAssociation(clusterVector, slidingFitResultMapPair, clusterAssociation))
            break;

        std::cout << "\033[31m" << "TrackMerge: Finding extrapolated hits..." << "\033[0m"  << std::endl;  
        ClusterToCaloHitListMap clusterToCaloHitListMap;
        this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, createdMainTrackClusters, clusterToCaloHitListMap);

        std::cout << "\033[31m" << "TrackMerge: Checking hits..." << "\033[0m"  << std::endl;     
        if (!this->AreExtrapolatedHitsGood(clusterAssociation, clusterToCaloHitListMap))
        {
            std::cout << "\033[31m" << "EXTRAPOLATED HITS ARE PANTS - ABORT" << "\033[0m"  << std::endl;
            this->ConsiderClusterAssociation(clusterAssociation, clusterVector, slidingFitResultMapPair);
            continue;
        }

        const ClusterList::const_iterator upstreamIter(std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), clusterAssociation.GetUpstreamCluster()));
        if (upstreamIter != createdMainTrackClusters.end())
            createdMainTrackClusters.erase(upstreamIter);

        const ClusterList::const_iterator downstreamIter(std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), clusterAssociation.GetDownstreamCluster()));
        if (downstreamIter != createdMainTrackClusters.end())
            createdMainTrackClusters.erase(downstreamIter);

        createdMainTrackClusters.push_back(this->CreateMainTrack(clusterAssociation, clusterToCaloHitListMap, pClusterList, clusterVector, slidingFitResultMapPair));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool TrackMergeRefinementAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair, ClusterPairAssociation &clusterAssociation) const
{
    bool foundAssociation(false);
    float maxLength(0.f);

    // ATTN: Find the associated cluster pair with the longest length sum 
    for (ClusterVector::const_iterator currentIter = clusterVector.begin(); currentIter != clusterVector.end(); ++currentIter)
    {
        const Cluster *const pCurrentCluster(*currentIter);

        const TwoDSlidingFitResultMap::const_iterator currentMicroFitIter(slidingFitResultMapPair.first->find(pCurrentCluster));
        if (currentMicroFitIter == slidingFitResultMapPair.first->end())
            return false;

        const TwoDSlidingFitResultMap::const_iterator currentMacroFitIter(slidingFitResultMapPair.second->find(pCurrentCluster));
        if (currentMacroFitIter == slidingFitResultMapPair.second->end())
            return false;

        for (ClusterVector::const_iterator testIter = std::next(currentIter); testIter != clusterVector.end(); ++testIter)
        {
            const Cluster *const pTestCluster(*testIter);

            const float lengthSum(LArClusterHelper::GetLength(pCurrentCluster) + LArClusterHelper::GetLength(pTestCluster));
            if ((lengthSum < maxLength) || (lengthSum < m_minClusterLengthSum))
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMicroFitIter(slidingFitResultMapPair.first->find(pTestCluster));
            if (testMicroFitIter == slidingFitResultMapPair.first->end())
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMacroFitIter(slidingFitResultMapPair.second->find(pTestCluster));
            if (testMacroFitIter == slidingFitResultMapPair.second->end())
                continue;

            const bool isCurrentUpstream(LArClusterHelper::SortByPosition(pCurrentCluster, pTestCluster));

            // ATTN: Ensure that clusters are not contained within one another
            const float currentMinLayerZ(currentMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()), currentMaxLayerZ(currentMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ());
            const float testMinLayerZ(testMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()), testMaxLayerZ(testMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ());

            if (((currentMinLayerZ > testMinLayerZ) && (currentMaxLayerZ < testMaxLayerZ)) || ((testMinLayerZ > currentMinLayerZ) && (testMaxLayerZ < currentMaxLayerZ)))
                continue;
            
            CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
            if (!this->GetClusterMergingCoordinates(currentMicroFitIter->second, currentMacroFitIter->second, testMacroFitIter->second, !isCurrentUpstream, currentMergePoint, currentMergeDirection) ||
                !this->GetClusterMergingCoordinates(testMicroFitIter->second, testMacroFitIter->second, currentMacroFitIter->second, isCurrentUpstream, testMergePoint, testMergeDirection))
            {
                continue;
            }

            // remember direction has changed
            if ((isCurrentUpstream && !this->AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection)) ||
                (!isCurrentUpstream && !this->AreClustersAssociated(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection)))
            {
                continue;
            }

            if (isCurrentUpstream)
            {
                clusterAssociation = ClusterPairAssociation(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection * (-1.f), pCurrentCluster, pTestCluster);
            }
            else
            {
                clusterAssociation = ClusterPairAssociation(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection * (-1.f), pTestCluster, pCurrentCluster);
            }

            foundAssociation = true;
            maxLength = lengthSum;
        }
    }

    if (foundAssociation)
    {
        ClusterList upCluster({clusterAssociation.GetUpstreamCluster()}), downCluster({clusterAssociation.GetDownstreamCluster()});
        const CartesianVector &upPoint(clusterAssociation.GetUpstreamMergePoint()), &downPoint(clusterAssociation.GetDownstreamMergePoint());
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &upCluster, "UPSTREAM CLUSTER", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &downCluster, "DOWNSTREAM CLUSTER", BLUE);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upPoint, "UP", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downPoint, "DOWN", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    else
    {
        std::cout << "DIDN'T FIND AN ASSOCIATION" << std::endl;
    }
    
    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackMergeRefinementAlgorithm::AreClustersAssociated(const CartesianVector &upstreamPoint, const CartesianVector &upstreamDirection, const CartesianVector &downstreamPoint,
    const CartesianVector &downstreamDirection) const
{
    if (downstreamPoint.GetZ() < upstreamPoint.GetZ())
        return false;

    const float separationDistance(std::sqrt(upstreamPoint.GetDistanceSquared(downstreamPoint)));
    if (separationDistance < m_minSeparationDistance)
        return false;

    if (upstreamDirection.GetCosOpeningAngle(downstreamDirection) < m_minDirectionDeviationCosAngle)
        return false;

    const CartesianVector predictedDownstreamPoint(upstreamPoint + (upstreamDirection * separationDistance));
    const float predictedDownstreamOffset((predictedDownstreamPoint - downstreamPoint).GetMagnitude());

    if (predictedDownstreamOffset > m_maxPredictedMergePointOffset)
        return false;

    const CartesianVector predictedUpstreamPoint(downstreamPoint - (downstreamDirection * separationDistance));
    const float predictedUpstreamOffset((predictedUpstreamPoint - upstreamPoint).GetMagnitude());

    if (predictedUpstreamOffset > m_maxPredictedMergePointOffset)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackMergeRefinementAlgorithm::GetExtrapolatedCaloHits(ClusterPairAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, const pandora::ClusterList &createdMainTrackClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    const float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));

    for (const Cluster *const pCluster : *pClusterList)
    {
        // ISOBEL THINK ABOUT THIS
        if ((std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), pCluster) != createdMainTrackClusters.end()) &&
            (pCluster != clusterAssociation.GetUpstreamCluster()) && (pCluster != clusterAssociation.GetDownstreamCluster()))
            continue;
        
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                CartesianVector closestPointToLine(LArHitWidthHelper::GetClosestPointToLine2D(upstreamPoint, clusterAssociation.GetConnectingLineDirection(), pCaloHit));

                //CHECKS IN X AND Z
                if ((closestPointToLine.GetX() < minX) || (closestPointToLine.GetX() > maxX) ||
                    (closestPointToLine.GetZ() < upstreamPoint.GetZ()) || (closestPointToLine.GetZ() > downstreamPoint.GetZ()))
                {
                    continue;
                }

                const float transverseDistanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(closestPointToLine - upstreamPoint).GetMagnitude());
                if (transverseDistanceFromLine > m_distanceFromLine)
                    continue;

                clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackMergeRefinementAlgorithm::AreExtrapolatedHitsGood(ClusterPairAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    CaloHitVector extrapolatedHitVector;
    for (const auto &entry : clusterToCaloHitListMap)
        extrapolatedHitVector.insert(extrapolatedHitVector.begin(), entry.second.begin(), entry.second.end());

    // ATTN: SORTED FROM UPSTREAM -> DOWNSTREAM POINT
    std::sort(extrapolatedHitVector.begin(), extrapolatedHitVector.end(), SortByDistanceAlongLine(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), m_hitWidthMode));

    //////////////////////////////////
    unsigned int count(0);
    for (const CaloHit *const pCaloHit : extrapolatedHitVector)
    {
        CartesianVector closestPoint(LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), pCaloHit));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint, std::to_string(count), YELLOW, 2);
        ++count;
    }
    const CartesianVector start(clusterAssociation.GetUpstreamMergePoint());
    const CartesianVector end(start + (clusterAssociation.GetConnectingLineDirection() * 150));
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "LINE", GREEN, 2, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    //////////////////////////////////    
    
    if (!this->AreExtrapolatedHitsNearMergePoints(extrapolatedHitVector, m_boundaryTolerance, clusterAssociation))
        return false;
    
    if (clusterToCaloHitListMap.empty())
        return true;

    if (!this->IsTrackContinuous(clusterAssociation, extrapolatedHitVector))
    {
        std::cout << "GAP IN HIT VECTOR" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    

    return true;
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackMergeRefinementAlgorithm::AreExtrapolatedHitsNearMergePoints(const CaloHitVector &extrapolatedHitVector, const float boundaryTolerance, 
    ClusterPairAssociation &clusterAssociation) const
{    
    if (extrapolatedHitVector.empty())
    {
        const float separationDistance((clusterAssociation.GetUpstreamMergePoint() - clusterAssociation.GetDownstreamMergePoint()).GetMagnitude());
        
        if (separationDistance > boundaryTolerance)
        {
            std::cout << "MERGE POINT TOO FAR AWAY FROM BOUNDARY & NO HITS" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            return false;
        }
        else
        {
            std::cout << "MERGE POINT CLOSE TO BOUNDARY & NO HITS" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            return true;
        }
    }

    const CaloHit *const closestCaloHit(extrapolatedHitVector.front()), *const furthestCaloHit(extrapolatedHitVector.back());
    
    const float distanceToUpstreamMergePoint(LArHitWidthHelper::GetClosestDistanceToPoint2D(closestCaloHit, clusterAssociation.GetUpstreamMergePoint()));
    const float distanceToDownstreamMergePoint(LArHitWidthHelper::GetClosestDistanceToPoint2D(furthestCaloHit, clusterAssociation.GetDownstreamMergePoint()));
    
    ////////////////////
    
    ClusterList upstream({clusterAssociation.GetUpstreamCluster()}), downstream({clusterAssociation.GetDownstreamCluster()});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &upstream, "UPSTREAM", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &downstream, "DOWNSTREAM", BLUE);
    CartesianVector closestPoint(LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), closestCaloHit));
    CartesianVector furthestPoint(LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), furthestCaloHit));
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint, "CLOSEST POINT", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &furthestPoint, "FURTHEST POINT", BLACK, 2);
    
    std::cout << "distanceToUpstreamMergePoint: " << distanceToUpstreamMergePoint << std::endl;
    std::cout << "distanceToDownstreamMergePoint: " << distanceToDownstreamMergePoint << std::endl;
    
    ////////////////////    
    
    if ((distanceToUpstreamMergePoint > boundaryTolerance) || (distanceToDownstreamMergePoint > boundaryTolerance))
    {
        std::cout << "failed cuts" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    return true;
}



//------------------------------------------------------------------------------------------------------------------------------------------

void TrackMergeRefinementAlgorithm::ConsiderClusterAssociation(const ClusterPairAssociation &clusterAssociation, pandora::ClusterVector &clusterVector,
    SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    const Cluster *const upstreamCluster(clusterAssociation.GetUpstreamCluster()), *const downstreamCluster(clusterAssociation.GetDownstreamCluster());
    const Cluster *const pConsideredCluster(upstreamCluster->GetNCaloHits() > downstreamCluster->GetNCaloHits() ? downstreamCluster : upstreamCluster);
    RemoveClusterFromContainers(pConsideredCluster, clusterVector, slidingFitResultMapPair);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackMergeRefinementAlgorithm::CreateMainTrack(const ClusterPairAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    ////////////////
    ClusterList upstreamTrack({clusterAssociation.GetUpstreamCluster()}), downstreamTrack({clusterAssociation.GetDownstreamCluster()});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &upstreamTrack, "UPSTREAM TRACK", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &downstreamTrack, "DOWNSTREAM TRACK", BLUE);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());    
    ///////////////
    
    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (const auto &entry : clusterToCaloHitListMap)
    {
        if ((entry.first != clusterAssociation.GetUpstreamCluster()) && (entry.first != clusterAssociation.GetDownstreamCluster()))
            showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ClusterList remnantClusterList;

    const Cluster *const pMainTrackCluster(this->RemoveOffAxisHitsFromTrack(clusterAssociation.GetUpstreamCluster(), clusterAssociation.GetUpstreamMergePoint(),
        false, clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second));
    const Cluster *const pClusterToDelete(this->RemoveOffAxisHitsFromTrack(clusterAssociation.GetDownstreamCluster(), clusterAssociation.GetDownstreamMergePoint(),
        true, clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pClusterToDelete));    

    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        this->AddHitsToMainTrack(pMainTrackCluster, pShowerCluster, caloHitsToMerge, clusterAssociation, remnantClusterList);
    }  

    ClusterList createdClusters;
    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList, createdClusters);

    ////////////////
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &createdClusters, "CREATED CLUSTERS", RED);
    ClusterList extendedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &extendedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());     
    ////////////////    

    // ATTN: Cleanup
    ClusterList modifiedClusters(showerClustersToFragment.begin(), showerClustersToFragment.end());
    modifiedClusters.push_back(clusterAssociation.GetUpstreamCluster()); modifiedClusters.push_back(clusterAssociation.GetDownstreamCluster());
    createdClusters.push_back(pMainTrackCluster);
    this->UpdateContainers(createdClusters, modifiedClusters, LArClusterHelper::SortByNHits, clusterVector, slidingFitResultMapPair);

    return pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackMergeRefinementAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLengthSum", m_minClusterLengthSum));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparationDistance", m_minSeparationDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPredictedMergePointOffset", m_maxPredictedMergePointOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFromLine", m_distanceFromLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BoundaryTolerance", m_boundaryTolerance));      

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------     

}// namespace lar_content
