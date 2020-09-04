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

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

using namespace pandora;

namespace lar_content
{

TrackExtensionRefinementAlgorithm::TrackExtensionRefinementAlgorithm() :
    m_maxLoopIterations(10),
    m_maxTrackDistanceToShowerBranch(5.f),    
    m_growingFitInitialLength(10.f),
    m_growingFitSegmentLength(5.0f),
    m_distanceToLine(1.0f),
    m_tpcBoundaryTolerance(2.f),
    m_mergePointBoundaryTolerance(2.f)
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

    // ATTN: Keep track of created main track clusters so their hits can be protected in future iterations
    ClusterList createdMainTrackClusters;
    for (bool isHigherXBoundary : { false, true })
    {
        m_isHigherXBoundary = isHigherXBoundary;

        /*
        const float nearestTPCBoundaryX(m_isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
        if ((std::fabs(nearestTPCBoundaryX - m_detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
            (std::fabs(nearestTPCBoundaryX - m_detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
        {
            continue;
        }
        */

        // ATTN: Keep track of clusters removed from clusterVector so that they can be added back in when considering the other endpoint
        ClusterList consideredClusters;
        unsigned int loopIterations(0);
        while(loopIterations < m_maxLoopIterations) 
        {
            ++loopIterations;

            ClusterEndpointAssociation clusterAssociation;
            if (!this->FindBestClusterAssociation(clusterVector, slidingFitResultMapPair, pClusterList, clusterAssociation))
                break;

            ClusterToCaloHitListMap clusterToCaloHitListMap;
            this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, createdMainTrackClusters, clusterToCaloHitListMap);
            
            if(!this->AreExtrapolatedHitsGood(clusterToCaloHitListMap, clusterAssociation))
            {
                this->ConsiderClusterAssociation(clusterAssociation.GetMainTrackCluster(), clusterAssociation.GetMainTrackCluster(), clusterVector, consideredClusters, slidingFitResultMapPair);
                continue;
            }

            const ClusterList::const_iterator iter(std::find(createdMainTrackClusters.begin(), createdMainTrackClusters.end(), clusterAssociation.GetMainTrackCluster()));
            if (iter != createdMainTrackClusters.end())
                createdMainTrackClusters.erase(iter);

            createdMainTrackClusters.push_back(this->CreateMainTrack(clusterAssociation, clusterToCaloHitListMap, pClusterList, clusterVector, slidingFitResultMapPair,
                consideredClusters));
        }

        if (!m_isHigherXBoundary)
            this->InitialiseContainers(&consideredClusters, SortByDistanceToTPCBoundary(m_tpcMaxXEdge), clusterVector, slidingFitResultMapPair);
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

bool TrackExtensionRefinementAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    const ClusterList *const pClusterList, ClusterEndpointAssociation &clusterAssociation) const
{
    const float nearestTPCBoundaryX(m_isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    
    // ATTN: This assumes that clusterVector has been sorted furthest cluster to TPC boundary -> closest
    for (const Cluster *const pCurrentCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCurrentCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCurrentCluster));

        const bool isEndUpstream = (std::fabs(microSlidingFitResult.GetGlobalMinLayerPosition().GetX() - nearestTPCBoundaryX) <
                                    std::fabs(microSlidingFitResult.GetGlobalMaxLayerPosition().GetX() - nearestTPCBoundaryX));
            
        CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isEndUpstream, clusterMergePoint, clusterMergeDirection))
            continue;

        // Reject clusters that do not cross TPC boundary
        if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
            continue;

        // If pCurrent cluster is contained within another close cluster it is likely to be part of a shower
        if (this->IsContained(pCurrentCluster, pClusterList))
            continue;

        // Reject clusters that do not fit criteria
        if(!this->DoesPassCriteria(microSlidingFitResult, clusterMergeDirection, isEndUpstream, pClusterList,  clusterMergePoint))
            continue;

        // ATTN: Temporarily set the other merge point to define extrapolate hits search region
        const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
        const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
        const CartesianVector extrapolatedHitsEndpoint(nearestTPCBoundaryX, 0.f, predictedIntercept + (predictedGradient * nearestTPCBoundaryX));

        if (isEndUpstream ? clusterMergePoint.GetZ() < extrapolatedHitsEndpoint.GetZ() : clusterMergePoint.GetZ() > extrapolatedHitsEndpoint.GetZ())
            continue;
                
        clusterAssociation = isEndUpstream ?
            ClusterEndpointAssociation(extrapolatedHitsEndpoint, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCurrentCluster, true) :
            ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedHitsEndpoint, clusterMergeDirection * (-1.f), pCurrentCluster, false);

        return true;
    }

    return false;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::IsContained(const Cluster *const pCurrentCluster, const ClusterList *const pClusterList) const
{
    const unsigned int currentInnerPseudoLayer(pCurrentCluster->GetInnerPseudoLayer()), currentOuterPseudoLayer(pCurrentCluster->GetOuterPseudoLayer());

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pCurrentCluster)
            continue;
                
        const unsigned int clusterInnerPseudoLayer(pCluster->GetInnerPseudoLayer()), clusterOuterPseudoLayer(pCluster->GetOuterPseudoLayer());

        if ((currentInnerPseudoLayer > clusterInnerPseudoLayer) && (currentOuterPseudoLayer < clusterOuterPseudoLayer))
        {
            if (LArClusterHelper::GetClosestDistance(pCurrentCluster, pCluster) < m_maxTrackDistanceToShowerBranch)
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void TrackExtensionRefinementAlgorithm::GetExtrapolatedCaloHits(const ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList,
    const ClusterList &createdMainTrackClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{    
    // Identify hits in the ROI
    ClusterList unavailableProtectedClusters;
    this->GetUnavailableProtectedClusters(clusterAssociation, createdMainTrackClusters, unavailableProtectedClusters);
        
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());

    ClusterToCaloHitListMap hitsInRegion;
    this->GetHitsInBoundingBox(upstreamPoint, downstreamPoint, pClusterList, hitsInRegion, unavailableProtectedClusters);
    
    ClusterVector clustersInRegion;
    for (const ClusterToCaloHitListMap::value_type &entry : hitsInRegion)
        clustersInRegion.push_back(entry.first);

    std::sort(clustersInRegion.begin(), clustersInRegion.end(), LArClusterHelper::SortByNHits);

    // Construct initial fit
    CartesianVector extrapolatedStartPosition(clusterAssociation.IsEndUpstream() ? downstreamPoint : upstreamPoint);
    CartesianVector extrapolatedDirection(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
    const CartesianVector clusterSubsetBoundary(extrapolatedStartPosition + (extrapolatedDirection * (-1.f) * m_growingFitInitialLength));

    ClusterToCaloHitListMap subsetFitHits;
    const ClusterList mainTrackClusterList({clusterAssociation.GetMainTrackCluster()});
    this->GetHitsInBoundingBox(extrapolatedStartPosition, clusterSubsetBoundary, &mainTrackClusterList, subsetFitHits);

    CartesianPointVector runningFitPositionVector;
    for (const CaloHit *const pCaloHit : subsetFitHits.begin()->second)
        runningFitPositionVector.push_back(pCaloHit->GetPositionVector());

    // Collect extrapolated hits by performing a running fit
    unsigned int count(0);
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
                    
                    if (!this->IsInLineSegment(extrapolatedStartPosition, extrapolatedEndPosition, hitPosition))
                        continue;

                    if (!this->IsCloseToLine(hitPosition, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
                        continue;

                    ++hitsCollected;
                    
                    runningFitPositionVector.push_back(hitPosition);
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
                }
            }
        }
        catch (const StatusCodeException &)
        {
            return;
        }

        ++count;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackExtensionRefinementAlgorithm::GetUnavailableProtectedClusters(const ClusterEndpointAssociation &clusterAssociation, const ClusterList &createdMainTrackClusters,
    ClusterList &unavailableProtectedClusters) const
{
    for (const Cluster *const pMainTrackCluster : createdMainTrackClusters)
    {
        if (pMainTrackCluster != clusterAssociation.GetMainTrackCluster())
            unavailableProtectedClusters.push_back(pMainTrackCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackExtensionRefinementAlgorithm::AreExtrapolatedHitsNearBoundaries(const CaloHitVector &extrapolatedHitVector, ClusterAssociation &clusterAssociation) const
{
    const float nearestTPCBoundaryX(m_isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    const bool isEndUpstream(std::fabs(clusterAssociation.GetUpstreamMergePoint().GetX() - nearestTPCBoundaryX) <
        std::fabs(clusterAssociation.GetDownstreamMergePoint().GetX() - nearestTPCBoundaryX));
    
    const CartesianVector &clusterMergePoint(isEndUpstream ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());

    if (extrapolatedHitVector.empty())
    {
        const float distanceFromTPCBoundary(std::fabs(clusterMergePoint.GetX() - nearestTPCBoundaryX));
        return (distanceFromTPCBoundary > m_mergePointBoundaryTolerance ? false : true);
    }

    const CaloHit *const furthestCaloHit(isEndUpstream ? extrapolatedHitVector.front() : extrapolatedHitVector.back());

    if (!this->IsNearBoundary(furthestCaloHit, CartesianVector(nearestTPCBoundaryX, 0.f, furthestCaloHit->GetPositionVector().GetZ()), m_tpcBoundaryTolerance))
        return false;

    const CaloHit *const closestCaloHit(isEndUpstream ? extrapolatedHitVector.back() : extrapolatedHitVector.front());

    if (!this->IsNearBoundary(closestCaloHit, clusterMergePoint, m_mergePointBoundaryTolerance))
        return false;

    // Reset extrapolated cluster merge point to be the projection of the furthest extrapolated hit
    const CartesianVector hitPosition(m_hitWidthMode ?
        LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), furthestCaloHit) :
        furthestCaloHit->GetPositionVector());
    const CartesianVector displacementVector(hitPosition - clusterMergePoint);
    const float signFactor(isEndUpstream ? -1.f : 1.f);
    const CartesianVector extrapolatedPoint(clusterMergePoint +
        (clusterAssociation.GetConnectingLineDirection() * std::fabs(displacementVector.GetDotProduct(clusterAssociation.GetConnectingLineDirection())) * signFactor * 1.0001));
    
    isEndUpstream ? clusterAssociation.SetUpstreamMergePoint(extrapolatedPoint) : clusterAssociation.SetDownstreamMergePoint(extrapolatedPoint);
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackExtensionRefinementAlgorithm::ConsiderClusterAssociation(const Cluster *const pOldConsideredCluster, const Cluster *const pNewConsideredCluster,
    ClusterVector &clusterVector, ClusterList &consideredClusters, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    RemoveClusterFromContainers(pOldConsideredCluster, clusterVector, slidingFitResultMapPair);
    consideredClusters.push_back(pNewConsideredCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackExtensionRefinementAlgorithm::CreateMainTrack(const ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterList &consideredClusters) const
{
    const Cluster *pMainTrackCluster(clusterAssociation.GetMainTrackCluster());
    const CartesianVector &clusterMergePoint(clusterAssociation.IsEndUpstream() ?
        clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());

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

    // ATTN: Cleanup containers - choose to not add created clusters back into containers
    ClusterList modifiedClusters(showerClustersToFragment.begin(), showerClustersToFragment.end());
    this->ConsiderClusterAssociation(clusterAssociation.GetMainTrackCluster(), pMainTrackCluster, clusterVector, consideredClusters, slidingFitResultMapPair);
    createdClusters.clear();
    this->UpdateContainers(createdClusters, modifiedClusters, SortByDistanceToTPCBoundary(m_isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge), clusterVector, slidingFitResultMapPair);

    return pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackExtensionRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLoopIterations", m_maxLoopIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackDistanceToShowerBranch", m_maxTrackDistanceToShowerBranch));    
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GrowingFitInitialLength", m_growingFitInitialLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GrowingFitSegmentLength", m_growingFitSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceToLine", m_distanceToLine));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TPCBoundaryTolerance", m_tpcBoundaryTolerance));       

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergePointBoundaryTolerance", m_mergePointBoundaryTolerance));    
    
    return TrackRefinementBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content

