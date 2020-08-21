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
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"

using namespace pandora;

namespace lar_content
{

TrackExtensionRefinementAlgorithm::TrackExtensionRefinementAlgorithm() :
    m_growingFitInitialLength(10.f),
    m_growingFitSegmentLength(5.0f),
    m_furthestDistanceToLine(10.f),
    m_closestDistanceToLine(1.0f),
    m_maxLoopIterations(10),
    m_distanceToLine(0.5f),
    m_boundaryTolerance(2.f)   
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
StatusCode TrackExtensionRefinementAlgorithm::Run()
{
     const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    ClusterVector clusterVector;
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    SlidingFitResultMapPair slidingFitResultMapPair({&microSlidingFitResultMap, &macroSlidingFitResultMap});

    this->InitialiseGeometry();
    
    this->InitialiseContainers(pClusterList, SortByDistanceToTPCBoundary(m_tpcMinXEdge), clusterVector, slidingFitResultMapPair);

    ClusterList createdMainTrackClusters;
    for (unsigned int isHigherXBoundary = 0; isHigherXBoundary < 2; ++isHigherXBoundary)
    {
        const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
        if ((std::fabs(nearestTPCBoundaryX - m_detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
            (std::fabs(nearestTPCBoundaryX - m_detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
        {
            continue;
        }

        ClusterList consideredClusters;
        unsigned int loopIterations(0);
        while(loopIterations < m_maxLoopIterations) 
        {
            ++loopIterations;

            ClusterEndpointAssociation clusterAssociation;
            if (!this->FindBestClusterAssociation(clusterVector, slidingFitResultMapPair, pClusterList, isHigherXBoundary, clusterAssociation))
                break;

            ClusterToCaloHitListMap clusterToCaloHitListMap;
            this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, createdMainTrackClusters, clusterToCaloHitListMap);

            if(!this->AreExtrapolatedHitsGood(clusterAssociation, clusterToCaloHitListMap, isHigherXBoundary))
            {
                this->ConsiderClusterAssociation(clusterAssociation.GetMainTrackCluster(), clusterAssociation.GetMainTrackCluster(), clusterVector, consideredClusters, slidingFitResultMapPair);
                continue;
            }
            
            std::cout << "\033[31m" << "Creating main track..." << "\033[0m"  << std::endl;
            const ClusterList::const_iterator iter(std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), clusterAssociation.GetMainTrackCluster()));
            if (iter != createdMainTrackClusters.end())
            {
                //std::cout << "GETTING CLUSTER OUT OF CREATED CLUSTERS" << std::endl;
                createdMainTrackClusters.erase(iter);
            }
            createdMainTrackClusters.push_back(this->CreateMainTrack(clusterAssociation, clusterToCaloHitListMap, pClusterList, clusterVector, slidingFitResultMapPair,
                consideredClusters, isHigherXBoundary));
        }

        if (isHigherXBoundary == 0)
            InitialiseContainers(&consideredClusters, SortByDistanceToTPCBoundary(m_tpcMaxXEdge), clusterVector, slidingFitResultMapPair);

        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &createdMainTrackClusters, "createdMainTrackClusters", GREEN);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }


    return STATUS_CODE_SUCCESS;
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
    
    const pandora::LArTPC *pLArTPC = &this->GetPandora().GetGeometry()->GetLArTPC();
    const LArTPCMap &larTPCMap(pPrimaryPandoraInstance->GetGeometry()->GetLArTPCMap());
    
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pSubLArTPC(mapEntry.second);
        m_detectorMinXEdge = std::min(m_detectorMinXEdge, pSubLArTPC->GetCenterX() - 0.5f * pSubLArTPC->GetWidthX());
        m_detectorMaxXEdge = std::max(m_detectorMaxXEdge, pSubLArTPC->GetCenterX() + 0.5f * pSubLArTPC->GetWidthX());

        //ATTN: Child & parent pandora instance TPCs have different addresses
        if (std::fabs(pSubLArTPC->GetCenterX() - pLArTPC->GetCenterX()) < std::numeric_limits<float>::epsilon())
            pLArTPC = pSubLArTPC;
    }

    m_tpcMinXEdge = pLArTPC->GetCenterX() - (pLArTPC->GetWidthX() * 0.5f);
    m_tpcMaxXEdge = pLArTPC->GetCenterX() + (pLArTPC->GetWidthX() * 0.5f);

    float &cpaXBoundary((pLArTPC->IsDriftInPositiveX() ? m_tpcMinXEdge : m_tpcMaxXEdge));

    if ((std::fabs(cpaXBoundary - m_detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
        (std::fabs(cpaXBoundary - m_detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
    {
        return;
    }
    
    const LArTPC *const pNeighboughTPC(&LArStitchingHelper::FindClosestTPC(*pPrimaryPandoraInstance, *pLArTPC, !pLArTPC->IsDriftInPositiveX()));
    const float gapSizeX(std::fabs(pNeighboughTPC->GetCenterX() - pLArTPC->GetCenterX()) - (pNeighboughTPC->GetWidthX() * 0.5f) - (pLArTPC->GetWidthX() * 0.5f));

    cpaXBoundary += gapSizeX * (pLArTPC->IsDriftInPositiveX() ? -0.5f : 0.5f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackExtensionRefinementAlgorithm::GetExtrapolatedCaloHits(ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList, const ClusterList &createdMainTrackClusters,
    ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    // Identify hits in the ROI
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
                CartesianVector hitPosition(m_hitWidthMode ?
                    LArHitWidthHelper::GetClosestPointToLine2D(upstreamPoint, clusterAssociation.GetConnectingLineDirection(), pCaloHit) :
                    pCaloHit->GetPositionVector());

                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;
                
                hitsInRegion[pCluster].push_back(pCaloHit);
            }
        }
    }

    ClusterVector clustersInRegion;
    for (const ClusterToCaloHitListMap::value_type &entry : hitsInRegion)
        clustersInRegion.push_back(entry.first);

    std::sort(clustersInRegion.begin(), clustersInRegion.end(), LArClusterHelper::SortByNHits);

    // Construct initial fit
    CartesianVector extrapolatedStartPosition(clusterAssociation.IsEndUpstream() ? downstreamPoint : upstreamPoint);
    CartesianVector extrapolatedDirection(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
    const CartesianVector clusterSubsetBoundary(extrapolatedStartPosition + (extrapolatedDirection * (-1.f) * m_growingFitInitialLength));

    minX = std::min(extrapolatedStartPosition.GetX(), clusterSubsetBoundary.GetX()); maxX = std::max(extrapolatedStartPosition.GetX(), clusterSubsetBoundary.GetX());
    const float minZ(std::min(extrapolatedStartPosition.GetZ(), clusterSubsetBoundary.GetZ())), maxZ(std::max(extrapolatedStartPosition.GetZ(), clusterSubsetBoundary.GetZ()));

    CartesianPointVector runningFitPositionVector;
    const OrderedCaloHitList &orderedCaloHitList(clusterAssociation.GetMainTrackCluster()->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < minZ) || (hitPosition.GetZ() > maxZ))
                continue;

            runningFitPositionVector.push_back(hitPosition);    
        }
    }

    // Collect extrapolated hits by performing a running fit
    unsigned int count(0);
    unsigned int segmentsWithoutHits(0);
    CartesianVector extrapolatedEndPosition(0.f, 0.f, 0.f);
    unsigned int hitsCollected(std::numeric_limits<int>::max());
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    while (hitsCollected)
    {
        hitsCollected = 0;

        try
        {
            const TwoDSlidingFitResult extrapolatedFit(&runningFitPositionVector, m_macroSlidingFitWindow, slidingFitPitch);
            
            if (count > 0)
            {
                extrapolatedStartPosition = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerPosition() : extrapolatedFit.GetGlobalMaxLayerPosition();
                extrapolatedDirection = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f) : extrapolatedFit.GetGlobalMaxLayerDirection();
            }
            
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);
            
            ////////////////
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedStartPosition, "start", RED, 2);
            ////////////////

            for (const Cluster *const pCluster : clustersInRegion)
            {
                for (const CaloHit *const pCaloHit : hitsInRegion.at(pCluster))
                {
                    // ATTN: To avoid counting same hit twice
                    const ClusterToCaloHitListMap::iterator iter(clusterToCaloHitListMap.find(pCluster));
                    if (iter != clusterToCaloHitListMap.end())
                    {
                        if (std::find(iter->second.begin(), iter->second.end(), pCaloHit) != iter->second.end())
                            continue;
                    }

                    CartesianVector hitPosition(m_hitWidthMode ?
                        LArHitWidthHelper::GetClosestPointToLine2D(extrapolatedStartPosition, extrapolatedDirection, pCaloHit) :
                        pCaloHit->GetPositionVector());
                    
                    if (!IsInLineSegment(extrapolatedStartPosition, extrapolatedEndPosition, hitPosition))
                        continue;

                    const float transverseDistanceFromLine(extrapolatedDirection.GetCrossProduct(hitPosition - extrapolatedStartPosition).GetMagnitude());
                    
                    if (transverseDistanceFromLine > m_distanceToLine)
                        continue;

                    ++hitsCollected;
                    
                    runningFitPositionVector.push_back(hitPosition);
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPointToLine, std::to_string(transverseDistanceFromLine), GREEN, 2);
                }
            }

            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedEndPosition, "end", RED, 2);
        }
        catch (const StatusCodeException &)
        {
            return;
        }

        ++count;

        if (!hitsCollected)
            ++segmentsWithoutHits;
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackExtensionRefinementAlgorithm::CreateMainTrack(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterList &consideredClusters, const bool isHigherXBoundary) const
{
    const Cluster *pMainTrackCluster(clusterAssociation.GetMainTrackCluster());
    const CartesianVector &clusterMergePoint(clusterAssociation.IsEndUpstream() ?
        clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());

    ////////////////
    /*
    ClusterList originalTrack({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &originalTrack, "ORIGINAL TRACK", BLACK);
    */
    ///////////////

    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (auto &entry : clusterToCaloHitListMap)
    {
        if (entry.first == pMainTrackCluster)
            continue;        

        const ClusterList::const_iterator iter(std::find(consideredClusters.begin(), consideredClusters.end(), entry.first));
        if (iter != consideredClusters.end())
            consideredClusters.erase(iter);
        
        showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ClusterList remnantClusterList;
    pMainTrackCluster = RemoveOffAxisHitsFromTrack(pMainTrackCluster, clusterMergePoint, clusterAssociation.IsEndUpstream(),
        clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second);
   
    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        this->AddHitsToMainTrack(pMainTrackCluster, pShowerCluster, caloHitsToMerge, clusterAssociation, remnantClusterList);
    }

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

    // ATTN: Cleanup containers - choose to not add created clusters back into containers
    ClusterList modifiedClusters(showerClustersToFragment.begin(), showerClustersToFragment.end());
    this->ConsiderClusterAssociation(clusterAssociation.GetMainTrackCluster(), pMainTrackCluster, clusterVector, consideredClusters, slidingFitResultMapPair);
    createdClusters.clear();
    this->UpdateContainers(createdClusters, modifiedClusters, SortByDistanceToTPCBoundary(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge), clusterVector, slidingFitResultMapPair);

    return pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackExtensionRefinementAlgorithm::ConsiderClusterAssociation(const Cluster *const pOldConsideredCluster, const Cluster *const pNewConsideredCluster, ClusterVector &clusterVector, ClusterList &consideredClusters, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    RemoveClusterFromContainers(pOldConsideredCluster, clusterVector, slidingFitResultMapPair);
    consideredClusters.push_back(pNewConsideredCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::AreExtrapolatedHitsGood(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const bool isHigherXBoundary) const
{
    CaloHitVector extrapolatedHitVector;
    for (const auto &entry : clusterToCaloHitListMap)
        extrapolatedHitVector.insert(extrapolatedHitVector.begin(), entry.second.begin(), entry.second.end());

    // ATTN: Extrapolated hit checks require extrapolatedHitVector to be ordered from upstream -> downstream merge point  
    std::sort(extrapolatedHitVector.begin(), extrapolatedHitVector.end(), SortByDistanceAlongLine(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), m_hitWidthMode));

    ////////////////////////////////////
    /*
    unsigned int count(0);
    for (const CaloHit *const pCaloHit : extrapolatedHitVector)
    {
        CartesianVector closestPoint(LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), pCaloHit));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint, std::to_string(count), GREEN, 2);
        ++count;
    }
    const CartesianVector start(clusterAssociation.GetUpstreamMergePoint());
    const CartesianVector end(start + (clusterAssociation.GetConnectingLineDirection() * 150));
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "LINE", GREEN, 2, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ////////////////////////////////////    
    
    if (!this->IsExtrapolatedEndpointNearBoundary(extrapolatedHitVector, isHigherXBoundary, clusterAssociation))
        return false;

    if (clusterToCaloHitListMap.empty())
        return true;

    // ISOBEL - remove
    if (!this->IsTrackContinuous(clusterAssociation, extrapolatedHitVector))
    {
        //std::cout << "GAP IN HIT VECTOR" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::IsExtrapolatedEndpointNearBoundary(const CaloHitVector &extrapolatedHitVector, const bool isHigherXBoundary, 
    ClusterEndpointAssociation &clusterAssociation) const
{
    const CartesianVector clusterMergePoint(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);

    ////////////////////
    /*
    ClusterList theCluster({clusterAssociation.GetMainTrackCluster()});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLACK);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "MERGE POINT", BLACK, 2);
    */
    ////////////////////  
    
    if (extrapolatedHitVector.empty())
    {
        const float distanceFromTPCBoundary(std::fabs(clusterMergePoint.GetX() - nearestTPCBoundaryX));
        
        if (distanceFromTPCBoundary > m_boundaryTolerance)
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

    const CaloHit *const furthestCaloHit(clusterAssociation.IsEndUpstream() ? extrapolatedHitVector.front() : extrapolatedHitVector.back());

    const CartesianVector furthestHitPosition(m_hitWidthMode ?
        LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), furthestCaloHit) :
        furthestCaloHit->GetPositionVector());

    const float distanceFromTPCBoundary(std::fabs(furthestHitPosition.GetX() - nearestTPCBoundaryX));

    if ((distanceFromTPCBoundary > m_boundaryTolerance))
    {
        //std::cout << "failed cuts" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    const CaloHit *const closestCaloHit(clusterAssociation.IsEndUpstream() ? extrapolatedHitVector.back() : extrapolatedHitVector.front());
    
    const float distanceToClusterMergePoint(m_hitWidthMode ?
        LArHitWidthHelper::GetClosestDistanceToPoint2D(closestCaloHit, clusterMergePoint) :
        (closestCaloHit->GetPositionVector() - clusterMergePoint).GetMagnitude());

    if (distanceToClusterMergePoint > m_boundaryTolerance)
    {
        //std::cout << "failed cuts" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    ////////////////////
    /*
    const CartesianVector &closestPoint(closestCaloHit->GetPositionVector());
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint, "CLOSEST POINT", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &furthestPoint, "FURTHEST POINT", RED, 2);
    
    std::cout << "distanceFromTPCBoundary: " << distanceFromTPCBoundary << std::endl;
    std::cout << "distance between merge and closest" << distanceToClusterMergePoint << std::endl;
    */
    ////////////////////    

    // Reset other cluster merge point to be the projection of the furthest extrapolated hit <--- ISOBEL: this doesn't need to be done? 
    const CartesianVector displacementVector(furthestHitPosition - clusterMergePoint);
    const float signFactor(clusterAssociation.IsEndUpstream() ? -1.f : 1.f);
    const CartesianVector extrapolatedPoint(clusterMergePoint +
        (clusterAssociation.GetConnectingLineDirection() * std::fabs(displacementVector.GetDotProduct(clusterAssociation.GetConnectingLineDirection())) * signFactor * 1.0001));
    
    clusterAssociation.IsEndUpstream() ? clusterAssociation.SetUpstreamMergePoint(extrapolatedPoint) : clusterAssociation.SetDownstreamMergePoint(extrapolatedPoint);

    /*
    const CartesianVector &upstream(clusterAssociation.GetUpstreamMergePoint());
    const CartesianVector &downstream(clusterAssociation.GetDownstreamMergePoint());
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstream, "UP", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstream, "DOWN", BLACK, 2);        
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
    return true;
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLoopIterations", m_maxLoopIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceToLine", m_distanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BoundaryTolerance", m_boundaryTolerance));       

    return TrackRefinementBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
