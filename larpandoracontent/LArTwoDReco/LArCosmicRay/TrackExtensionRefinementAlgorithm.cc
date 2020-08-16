/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackExtensionRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the track refinement class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackExtensionRefinementAlgorithm.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

TrackExtensionRefinementAlgorithm::TrackExtensionRefinementAlgorithm() :
    m_growingFitInitialLength(20.f),
    m_growingFitSegmentLength(5.0f),
    m_furthestDistanceToLine(10.f),
    m_closestDistanceToLine(0.5f)        
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
StatusCode TrackExtensionRefinementAlgorithm::Run()
{
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    //std::cout << "IF TRACK IN EM SHOWER REMEMBER YOU CHANGED THE DIRECTION!" << std::endl;
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    ClusterVector clusterVector;
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    SlidingFitResultMapPair slidingFitResultMapPair({&microSlidingFitResultMap, &macroSlidingFitResultMap});

    this->InitialiseGeometry();
    
    this->InitialiseContainers(pClusterList, clusterVector, slidingFitResultMapPair);

    /*
    if (std::fabs(m_tpcMinXEdge - 3.24704) > 0.1)
        return STATUS_CODE_SUCCESS;
    */
    ClusterList createdMainTrackClusters;
    for (unsigned int isHigherXBoundary = 0; isHigherXBoundary < 2; ++isHigherXBoundary)
    {

        //std::cout << "\033[31m" << "isHigherXBoundary: " << isHigherXBoundary << "\033[0m"  << std::endl;
        
        const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
        if ((std::fabs(nearestTPCBoundaryX - m_detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
            (std::fabs(nearestTPCBoundaryX - m_detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
        {
            continue;
        }
        
        unsigned int loopIterations(0);
        ClusterList consideredClusters;
        while(loopIterations < 10) 
        {
            ++loopIterations;

            std::sort(clusterVector.begin(), clusterVector.end(), SortByDistanceToTPCBoundary(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge));

            //std::cout << "\033[31m" << "Finding best cluster association..." << "\033[0m"  << std::endl;
            ClusterEndpointAssociation clusterAssociation;
            if (!this->FindBestClusterAssociation(clusterVector, slidingFitResultMapPair, clusterAssociation, pClusterList, isHigherXBoundary))
                break;

            //std::cout << "\033[31m" << "Finding extrapolated hits..." << "\033[0m"  << std::endl;            
            ClusterToCaloHitListMap clusterToCaloHitListMap;
            this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, createdMainTrackClusters, clusterToCaloHitListMap);

            //std::cout << "\033[31m" << "Asking whether hits make sense..." << "\033[0m"  << std::endl;              
            if(!this->AreExtrapolatedHitsGood(clusterAssociation, clusterToCaloHitListMap, isHigherXBoundary))
            {
                //std::cout << "\033[31m" << "EXTRAPOLATED HITS ARE PANTS - ABORT" << "\033[0m"  << std::endl;
                this->ConsiderClusterAssociation(clusterAssociation, clusterVector, consideredClusters, slidingFitResultMapPair);
                continue;
            }

            //std::cout << "\033[31m" << "Creating main track..." << "\033[0m"  << std::endl;
            const ClusterList::const_iterator iter(std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), clusterAssociation.GetMainTrackCluster()));
            if (iter != createdMainTrackClusters.end())
            {
                //std::cout << "GETTING CLUSTER OUT OF CREATED CLUSTERS" << std::endl;
                createdMainTrackClusters.erase(iter);
            }
            createdMainTrackClusters.push_back(this->CreateMainTrack(clusterAssociation, clusterToCaloHitListMap, pClusterList, clusterVector, slidingFitResultMapPair, consideredClusters));
        }

        //std::cout << "\033[31m" << "Adding considered clusters back into the cluster vector..." << "\033[0m"  << std::endl;          
        //std::cout << "considered clusters size: " << consideredClusters.size() << std::endl;
        if ((isHigherXBoundary == 0) && (!consideredClusters.empty()))
            InitialiseContainers(&consideredClusters, clusterVector, slidingFitResultMapPair);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &createdMainTrackClusters, "createdMainTrackClusters", GREEN);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------  

void TrackExtensionRefinementAlgorithm::GetExtrapolatedCaloHits(ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList, const ClusterList &createdMainTrackClusters,
    ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{

    /*
    for (const Cluster *const pCluster : createdMainTrackClusters)
    {
        if (pCluster == clusterAssociation.GetMainTrackCluster())
            continue;
        
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "CANNOT HAVE", BLACK, 2);
            }
        }
    }
    */

    
    // Look for clusters in the region of interest
    ClusterToCaloHitListMap hitsInRegion;
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));
    for (const Cluster *const pCluster : *pClusterList)
    {

        if ((std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), pCluster) != createdMainTrackClusters.end()) && (pCluster != clusterAssociation.GetMainTrackCluster()))
            continue;
        
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                /*
                if(!IsInLineSegment(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetDownstreamMergePoint(), hitPosition))
                    continue;

                */
                //CHECKS IN X AND Z
                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;               
                
                hitsInRegion[pCluster].push_back(pCaloHit);
            }
        }
    }

    // ATTN: hitsInRegion is ordered (it isnt)
    ClusterVector clustersInRegion;
    for (const ClusterToCaloHitListMap::value_type &entry : hitsInRegion)
        clustersInRegion.push_back(entry.first);

    std::sort(clustersInRegion.begin(), clustersInRegion.end(), LArClusterHelper::SortByNHits);

    
    CartesianVector extrapolatedStartPosition(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());
    CartesianVector extrapolatedDirection(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
    const CartesianVector clusterSubsetBoundary(extrapolatedStartPosition + (extrapolatedDirection * (-1.f) * m_growingFitInitialLength));

    minX = std::min(extrapolatedStartPosition.GetX(), clusterSubsetBoundary.GetX()); maxX = std::max(extrapolatedStartPosition.GetX(), clusterSubsetBoundary.GetX());
    const float minZ(std::min(extrapolatedStartPosition.GetZ(), clusterSubsetBoundary.GetZ())), maxZ(std::max(extrapolatedStartPosition.GetZ(), clusterSubsetBoundary.GetZ()));

    CartesianPointVector hitPositionVector;
    const OrderedCaloHitList &orderedCaloHitList(clusterAssociation.GetMainTrackCluster()->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < minZ) || (hitPosition.GetZ() > maxZ))
                continue;

            hitPositionVector.push_back(hitPosition);    
        }
    }
    
    unsigned int count(0);
    //unsigned int count(1);
    unsigned int hitsCollected(std::numeric_limits<int>::max());
    CartesianVector extrapolatedEndPosition(0.f, 0.f, 0.f);
    unsigned int segmentsWithoutHits(0);

    int slidingFitWindow(m_microSlidingFitWindow);
    while (segmentsWithoutHits < 2)
    {
        hitsCollected = 0;

        float hitWidthCount(0.f);
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const TwoDSlidingFitResult extrapolatedFit(&hitPositionVector, slidingFitWindow, slidingFitPitch); //was micro m_macroSlidingFitWindow 
            
            if (count > 0)
            {
                extrapolatedStartPosition = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerPosition() : extrapolatedFit.GetGlobalMaxLayerPosition();
                //extrapolatedStartPosition = extrapolatedEndPosition;
                extrapolatedDirection = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f) : extrapolatedFit.GetGlobalMaxLayerDirection();
            }
            
            
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);
            
            if (segmentsWithoutHits == 1)
                extrapolatedEndPosition += (extrapolatedDirection * m_growingFitSegmentLength);
            
            const float gradient((extrapolatedEndPosition.GetZ() - extrapolatedStartPosition.GetZ()) / (extrapolatedEndPosition.GetX() - extrapolatedStartPosition.GetX()));
                        
            ////////////////
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedStartPosition, "start", RED, 2);
            ////////////////

            float distanceToLine(std::fabs(slidingFitWindow - m_microSlidingFitWindow) < std::numeric_limits<float>::epsilon() ? m_closestDistanceToLine : 2.0);
            for (const Cluster *const pCluster : clustersInRegion)
            {
                const CaloHitList &theCaloHits(hitsInRegion.at(pCluster));
                for (const CaloHit *const pCaloHit : theCaloHits)
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

                    // ATTN: To avoid counting same hit twice
                    const ClusterToCaloHitListMap::iterator iter(clusterToCaloHitListMap.find(pCluster));
                    if (iter != clusterToCaloHitListMap.end())
                    {
                        if (std::find(iter->second.begin(), iter->second.end(), pCaloHit) != iter->second.end())
                            continue;
                    }

                    if (!IsInLineSegment(extrapolatedStartPosition, extrapolatedEndPosition, hitPosition))
                        continue;

                    const float &hitWidth(pCaloHit->GetCellSize1());

                    const CartesianVector hitHighEdge(hitPosition.GetX() + (hitWidth * 0.5f), 0, hitPosition.GetZ());
                    const CartesianVector hitLowEdge(hitPosition.GetX() - (hitWidth * 0.5f), 0, hitPosition.GetZ());
                    
                    const float highEdgeDistanceFromLine((extrapolatedEndPosition - extrapolatedStartPosition).GetCrossProduct(hitHighEdge - extrapolatedStartPosition).GetMagnitude());
                    const float lowEdgeDistanceFromLine((extrapolatedEndPosition - extrapolatedStartPosition).GetCrossProduct(hitLowEdge - extrapolatedStartPosition).GetMagnitude());

                    //if ((highEdgeDistanceFromLine < 5.f) || (lowEdgeDistanceFromLine < 5.f))
                    //{
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitHighEdge, std::to_string(highEdgeDistanceFromLine), VIOLET, 2);
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitLowEdge, std::to_string(lowEdgeDistanceFromLine), VIOLET, 2);
                        //}
                    
                    if ((highEdgeDistanceFromLine > m_furthestDistanceToLine) || (lowEdgeDistanceFromLine > m_furthestDistanceToLine))
                        continue;

                    float xOnLine(((hitPosition.GetZ() - extrapolatedStartPosition.GetZ()) / gradient) + extrapolatedStartPosition.GetX());

                    if (((hitHighEdge.GetX() > xOnLine) && (hitLowEdge.GetX() > xOnLine)) ||
                        ((hitHighEdge.GetX() < xOnLine) && (hitLowEdge.GetX() < xOnLine)))
                    {
                    
                        if (!((highEdgeDistanceFromLine < distanceToLine) || (lowEdgeDistanceFromLine < distanceToLine)))
                            continue;
                    }

                    hitWidthCount += hitWidth;
                    ++hitsCollected;
                    hitPositionVector.push_back(hitPosition);
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "added", GREEN, 2);
                }
            }

            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedEndPosition, "end", RED, 2);
        }
        catch (const StatusCodeException &)
        {
            std::cout << "cannot make fit" << std::endl;
            return;
        }

        ++count;

        if (!hitsCollected)
            ++segmentsWithoutHits;

        if (hitsCollected)
        {
            //std::cout << "sparseness: " << (hitWidthCount / hitsCollected) << std::endl;
            slidingFitWindow = (hitWidthCount / hitsCollected) > 0.6 ? m_macroSlidingFitWindow : m_microSlidingFitWindow;
        }
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackExtensionRefinementAlgorithm::CreateMainTrack(ClusterEndpointAssociation &clusterEndpointAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterList &consideredClusters) const
{
    const Cluster *pMainTrackCluster(clusterEndpointAssociation.GetMainTrackCluster());
    const CartesianVector &clusterMergePoint(clusterEndpointAssociation.IsEndUpstream() ?
        clusterEndpointAssociation.GetDownstreamMergePoint() : clusterEndpointAssociation.GetUpstreamMergePoint());

   
    ////////////////
    //ClusterList originalTrack({pMainTrackCluster});
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &originalTrack, "ORIGINAL TRACK", BLACK);
    ///////////////

    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (auto &entry : clusterToCaloHitListMap)
    {
        if (entry.first == pMainTrackCluster)
            continue;        

        const ClusterList::const_iterator iter(std::find(consideredClusters.begin(), consideredClusters.end(), entry.first));
        if (iter != consideredClusters.end())
        {
            std::cout << "ISOBEL MODIFYING A CLUSTER THAT WAS IN CONSIDERED CLUSTERS" << std::endl;
            consideredClusters.erase(iter);
        }
        
        showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ////////////////
    /*
    for (auto &entry : showerClustersToFragment)
    {
        ClusterList jam({entry});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &jam, "SHOWER CLUSTER", VIOLET);
    }
    */
    ////////////////

    ClusterList remnantClusterList;
    pMainTrackCluster = RemoveOffAxisHitsFromTrack(pMainTrackCluster, clusterMergePoint, clusterEndpointAssociation.IsEndUpstream(),
        clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second);
   
    ////////////////
    /*
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    ClusterList refinedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &refinedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT CLUSTERS", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ////////////////
   
    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        this->AddHitsToMainTrack(pMainTrackCluster, pShowerCluster, caloHitsToMerge, clusterEndpointAssociation, remnantClusterList);
    }

    ////////////////
    /*
    ClusterList extendedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &extendedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT CLUSTERS", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());    
    */
    ////////////////   

    ClusterList createdClusters;
    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList, createdClusters);

    ////////////////
    /*
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &createdClusters, "CREATED CLUSTERS", RED);
    ClusterList extendedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &extendedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());     
    */
    ////////////////

    //ATTN: Update references to pMainTrackCluster
    //IT IS REALLY IMPORTANT THAT THIS IS DONE FIRST I.E. DELETE NULL POINTERS FROM MAP BEFORE ANY ARE ADDED IN AS THEY CAN BECOME THE SAME
    
    // ATTN: Cleanup
    ClusterList modifiedClusters(showerClustersToFragment.begin(), showerClustersToFragment.end());
    this->ConsiderClusterAssociation(clusterEndpointAssociation, clusterVector, consideredClusters, slidingFitResultMapPair);
    consideredClusters.push_back(pMainTrackCluster);
    createdClusters.clear();
    this->UpdateContainers(createdClusters, modifiedClusters, clusterVector, slidingFitResultMapPair);

    return pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
void TrackExtensionRefinementAlgorithm::UpdateAfterMainTrackModification(const Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation,
    SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    const Cluster *const pDeletedTrackCluster(clusterEndpointAssociation.GetMainTrackCluster());
    
    // REPLACE IN MICRO FIT MAP
    const TwoDSlidingFitResultMap::const_iterator microFitToDeleteIter(slidingFitResultMapPair.first->find(pDeletedTrackCluster));
    if (microFitToDeleteIter == slidingFitResultMapPair.first->end())
    {
        std::cout << "ISOBEL THIS IS BAD" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    slidingFitResultMapPair.first->erase(microFitToDeleteIter);

    // REPLACE IN MACRO FIT MAP
    const TwoDSlidingFitResultMap::const_iterator macroFitToDeleteIter(slidingFitResultMapPair.second->find(pDeletedTrackCluster));
    if (macroFitToDeleteIter == slidingFitResultMapPair.second->end())
    {
        std::cout << "ISOBEL THIS IS BAD" << std::endl;
        throw STATUS_CODE_FAILURE;
    }

    slidingFitResultMapPair.second->erase(macroFitToDeleteIter);


}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

void TrackExtensionRefinementAlgorithm::ConsiderClusterAssociation(const ClusterEndpointAssociation &clusterAssociation, ClusterVector &clusterVector, ClusterList &/*consideredClusters*/,
    SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    const Cluster *const pConsideredCluster(clusterAssociation.GetMainTrackCluster());
    RemoveClusterFromContainers(pConsideredCluster, clusterVector, slidingFitResultMapPair);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::AreExtrapolatedHitsGood(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const bool isHigherXBoundary) const
{
    const float m_boundaryTolerance(2.f);

    CaloHitVector extrapolatedHitVector;
    for (const auto &entry : clusterToCaloHitListMap)
        extrapolatedHitVector.insert(extrapolatedHitVector.begin(), entry.second.begin(), entry.second.end());

    // ATTN: SORTED FROM UPSTREAM -> DOWNSTREAM POINT
    std::sort(extrapolatedHitVector.begin(), extrapolatedHitVector.end(), SortByDistanceAlongLine(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection()));
    
    ////////////////////
    /*
    ClusterList theJam({clusterAssociation.GetMainTrackCluster()});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theJam, "MAIN TRACK", BLUE);
    for (const CaloHit *const pCaloHit : extrapolatedHitVector)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "EXTRAP", GREEN, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ///////////////    
    
    if (!this->IsExtrapolatedEndpointNearBoundary(extrapolatedHitVector, isHigherXBoundary, m_boundaryTolerance, clusterAssociation))
        return false;

    if (clusterToCaloHitListMap.empty())
        return true;

    const float averageHitSeparation(this->GetAverageHitSeparation(extrapolatedHitVector));
    std::cout << "averageHitSeparation: " << averageHitSeparation << std::endl;
    
    if (!this->IsTrackContinuous(clusterAssociation, extrapolatedHitVector, m_maxTrackGaps, m_lineSegmentLength))
    {
        std::cout << "GAP IN HIT VECTOR" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //consideredClusters.push_back(clusterAssociation.GetMainTrackCluster());                 
        return false;
    }

    return true;

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::IsExtrapolatedEndpointNearBoundary(const CaloHitVector &extrapolatedHitVector, const bool isHigherXBoundary, const float boundaryTolerance, 
    ClusterEndpointAssociation &clusterAssociation) const
{
    const CartesianVector clusterMergePoint(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);

    ////////////////////    
    //ClusterList theCluster({clusterAssociation.GetMainTrackCluster()});
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLACK);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "MERGE POINT", BLACK, 2);
    ////////////////////  
    
    if (extrapolatedHitVector.empty())
    {
        const float distanceFromTPCBoundary(std::fabs(clusterMergePoint.GetX() - nearestTPCBoundaryX));
        
        if (distanceFromTPCBoundary > boundaryTolerance)
        {
            //std::cout << "MERGE POINT TOO FAR AWAY FROM BOUNDARY & NO HITS" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            return false;
        }
        else
        {
            //std::cout << "MERGE POINT CLOSE TO BOUNDARY & NO HITS" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            return true;
        }
    }

    const CartesianVector &closestPoint(clusterAssociation.IsEndUpstream() ? extrapolatedHitVector.back()->GetPositionVector() : extrapolatedHitVector.front()->GetPositionVector());
    const CartesianVector &furthestPoint(clusterAssociation.IsEndUpstream() ? extrapolatedHitVector.front()->GetPositionVector() : extrapolatedHitVector.back()->GetPositionVector());
    
    const float distanceFromTPCBoundary(std::fabs(furthestPoint.GetX() - nearestTPCBoundaryX));

    ////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint, "CLOSEST POINT", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &furthestPoint, "FURTHEST POINT", RED, 2);
    std::cout << "distanceFromTPCBoundary: " << distanceFromTPCBoundary << std::endl;
    std::cout << "distance between merge and closest" << (clusterMergePoint - closestPoint).GetMagnitude() << std::endl;
    */    
    ////////////////////    
    
    if ((distanceFromTPCBoundary > boundaryTolerance) || ((clusterMergePoint - closestPoint).GetMagnitude() > 2.f))
    {
        //std::cout << "failed cuts" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    clusterAssociation.IsEndUpstream() ? clusterAssociation.SetUpstreamMergePoint(furthestPoint) : clusterAssociation.SetDownstreamMergePoint(furthestPoint);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackExtensionRefinementAlgorithm::GetAverageHitSeparation(const CaloHitVector &orderedCaloHitVector) const
{
    const CaloHit *pPreviousCaloHit(orderedCaloHitVector.front());

    float separationSum(0.f);
    for (const CaloHit *const pCaloHit : orderedCaloHitVector)
    {
        if (pCaloHit == pPreviousCaloHit)
            continue;

        separationSum += std::sqrt(pCaloHit->GetPositionVector().GetDistanceSquared(pPreviousCaloHit->GetPositionVector()));
        pPreviousCaloHit = pCaloHit;
    }

    return (separationSum / orderedCaloHitVector.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackExtensionRefinementAlgorithm::InitialiseGeometry()
{
    const Pandora *pPrimaryPandoraInstance;
    try
    {
        pPrimaryPandoraInstance = MultiPandoraApi::GetPrimaryPandoraInstance(&this->GetPandora());
    }
    catch (const StatusCodeException &)
    {
        pPrimaryPandoraInstance = &this->GetPandora();
    }
           
    m_detectorMinXEdge = std::numeric_limits<float>::max();
    m_detectorMaxXEdge = -std::numeric_limits<float>::max();
    m_pLArTPC = &this->GetPandora().GetGeometry()->GetLArTPC();
    
    const LArTPCMap &larTPCMap(pPrimaryPandoraInstance->GetGeometry()->GetLArTPCMap());
    
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pSubLArTPC(mapEntry.second);
        m_detectorMinXEdge = std::min(m_detectorMinXEdge, pSubLArTPC->GetCenterX() - 0.5f * pSubLArTPC->GetWidthX());
        m_detectorMaxXEdge = std::max(m_detectorMaxXEdge, pSubLArTPC->GetCenterX() + 0.5f * pSubLArTPC->GetWidthX());

        //ATTN: Child & parent pandora instance TPCs have different addresses
        if (std::fabs(pSubLArTPC->GetCenterX() - m_pLArTPC->GetCenterX()) < std::numeric_limits<float>::epsilon())
            m_pLArTPC = pSubLArTPC;
    }

    m_tpcMinXEdge = m_pLArTPC->GetCenterX() - (m_pLArTPC->GetWidthX() * 0.5f);
    m_tpcMaxXEdge = m_pLArTPC->GetCenterX() + (m_pLArTPC->GetWidthX() * 0.5f);

    float &cpaXBoundary((m_pLArTPC->IsDriftInPositiveX() ? m_tpcMinXEdge : m_tpcMaxXEdge));

    if ((std::fabs(cpaXBoundary - m_detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
        (std::fabs(cpaXBoundary - m_detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
    {
        /*
        const CartesianVector lowXPoint(m_tpcMinXEdge, m_pLArTPC->GetCenterY(), m_pLArTPC->GetCenterZ());
        const CartesianVector highXPoint(m_tpcMaxXEdge, m_pLArTPC->GetCenterY(), m_pLArTPC->GetCenterZ());     
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lowXPoint, "lowXPoint", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &highXPoint, "highXPoint", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return;
    }
    
    const LArTPC *const pNeighboughTPC(&LArStitchingHelper::FindClosestTPC(*pPrimaryPandoraInstance, *m_pLArTPC, !m_pLArTPC->IsDriftInPositiveX()));
    const float gapSizeX(std::fabs(pNeighboughTPC->GetCenterX() - m_pLArTPC->GetCenterX()) - (pNeighboughTPC->GetWidthX() * 0.5f) - (m_pLArTPC->GetWidthX() * 0.5f));

    cpaXBoundary += gapSizeX * (m_pLArTPC->IsDriftInPositiveX() ? -0.5f : 0.5f);

    /*
    const CartesianVector lowXPoint(m_tpcMinXEdge, m_pLArTPC->GetCenterY(), m_pLArTPC->GetCenterZ());
    const CartesianVector highXPoint(m_tpcMaxXEdge, m_pLArTPC->GetCenterY(), m_pLArTPC->GetCenterZ());

    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lowXPoint, "lowXPoint", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &highXPoint, "highXPoint", BLUE, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
}




//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackExtensionRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GrowingFitInitialLength", m_growingFitInitialLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GrowingFitSegmentLength", m_growingFitSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FurthestDistanceToLine", m_furthestDistanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClosestDistanceToLine", m_closestDistanceToLine));    

    return TrackRefinementBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::SortByDistanceToTPCBoundary::operator() (const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int lhsInnerPseudoLayer(pLhs->GetInnerPseudoLayer()), lhsOuterPseudoLayer(pLhs->GetOuterPseudoLayer());
    const float lhsInnerX(pLhs->GetCentroid(lhsInnerPseudoLayer).GetX()), lhsOuterX(pLhs->GetCentroid(lhsOuterPseudoLayer).GetX());

    const unsigned int rhsInnerPseudoLayer(pRhs->GetInnerPseudoLayer()), rhsOuterPseudoLayer(pRhs->GetOuterPseudoLayer());
    const float rhsInnerX(pRhs->GetCentroid(rhsInnerPseudoLayer).GetX()), rhsOuterX(pRhs->GetCentroid(rhsOuterPseudoLayer).GetX());       

    const float lhsFurthestDistance(std::max(std::fabs(lhsInnerX - m_tpcXBoundary), std::fabs(lhsOuterX - m_tpcXBoundary)));
    const float rhsFurthestDistance(std::max(std::fabs(rhsInnerX - m_tpcXBoundary), std::fabs(rhsOuterX - m_tpcXBoundary)));

    // ATTN: Order from furthest away to closest
    return (lhsFurthestDistance > rhsFurthestDistance);
}


    
}// namespace lar_content



/*
            /////////////////////

            std::cout  << "\033[31m" << "After hits collected..." << "\033[0m" << std::endl;
            CaloHitVector extrapolatedCaloHitVector;
            for (const auto &entry : clusterToCaloHitListMap)
                extrapolatedCaloHitVector.insert(extrapolatedCaloHitVector.begin(), entry.second.begin(), entry.second.end());
            
            for (auto &entry : extrapolatedCaloHitVector)
            {
                const CartesianVector &hitPosition(entry->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "EXTRAPOLATED HIT", GREEN, 2);
            }

            //////////////////////////////
*/
