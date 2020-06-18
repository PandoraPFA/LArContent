/**
 *  @file   TrackInEMShowerAlgorithm.cc
 *
 *  @brief  Implementation of the track in em shower algorithm class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/TrackInEMShowerAlgorithm.h"

using namespace pandora;

namespace lar_content
{
    
TrackInEMShowerAlgorithm::TrackInEMShowerAlgorithm() :
    m_minCaloHits(25), 
    m_minSeparationDistance(27.f),
    m_maxPredictedMergePointOffset(5.f),
    m_slidingFitWindow(20),
    m_mergePointMinCosAngleDeviation(0.9995f),
    m_minClusterLengthSum(75.f),
    m_minDirectionDeviationCosAngle(0.99f),
    m_distanceFromLine(0.35f),
    m_maxTrackGaps(2),
    m_lineSegmentLength(7.f),
    m_maxHitDistanceFromCluster(4.f),
    m_maxHitSeparationForConnectedCluster(4.f),
    m_maxDistanceFromMainTrack(0.75f),
    m_maxMainLoopIterations(5),
    m_macroSlidingFitWindow(1000),
    m_stableRegionClusterFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackInEMShowerAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap microSlidingFitResultMap;
    TwoDSlidingFitResultMap macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    // ATTN: Find longest length association and conitnue merging until no further merges can be made
    bool mergeMade(false);
    unsigned int loopIterations(0);
    do
    {
        ++loopIterations;
        
        ClusterAssociation clusterAssociation;
        if(!this->FindBestClusterAssociation(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap, clusterAssociation))
            break;

        CaloHitVector extrapolatedCaloHitVector;
        ClusterToCaloHitListMap clusterToCaloHitListMap;
        this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, extrapolatedCaloHitVector, clusterToCaloHitListMap);
        
        if (extrapolatedCaloHitVector.empty())
            break;
        
        if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
            break;
        
        this->CreateMainTrack(clusterAssociation, extrapolatedCaloHitVector, clusterToCaloHitListMap, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector);

        mergeMade = true;
        
    }
    while(mergeMade && (loopIterations < m_maxMainLoopIterations));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector) const
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

void TrackInEMShowerAlgorithm::InitialiseSlidingFitResultMaps(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            TwoDSlidingFitResult microSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);
            TwoDSlidingFitResult macroSlidingFitResult(pCluster, m_macroSlidingFitWindow, slidingFitPitch);

            (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, microSlidingFitResult));
            (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, macroSlidingFitResult));
        }
        catch (StatusCodeException &) {}
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::FindBestClusterAssociation(ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
    const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation) const
{
    bool foundAssociation(false);
    float maxLength(0.f);

    // ATTN: Find the associated cluster pair with the longest length sum
    for (ClusterVector::const_iterator currentIter = clusterVector.begin(); currentIter != clusterVector.end(); ++currentIter)
    {
        const Cluster *const pCurrentCluster(*currentIter);
        
        const TwoDSlidingFitResultMap::const_iterator currentMicroFitIter(microSlidingFitResultMap.find(pCurrentCluster));
        if (currentMicroFitIter == microSlidingFitResultMap.end())
            return false;

        const TwoDSlidingFitResultMap::const_iterator currentMacroFitIter(macroSlidingFitResultMap.find(pCurrentCluster));
        if (currentMacroFitIter == macroSlidingFitResultMap.end())
            return false;
   
        for (ClusterVector::const_iterator testIter = std::next(currentIter); testIter != clusterVector.end(); ++testIter)
        {
            const Cluster *const pTestCluster(*testIter);

            const float lengthSum(LArClusterHelper::GetLength(pCurrentCluster) + LArClusterHelper::GetLength(pTestCluster)); 
            if ((lengthSum < maxLength) || (lengthSum < m_minClusterLengthSum))
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMicroFitIter(microSlidingFitResultMap.find(pTestCluster));
            if (testMicroFitIter == microSlidingFitResultMap.end())
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMacroFitIter(macroSlidingFitResultMap.find(pTestCluster));
            if (testMacroFitIter == macroSlidingFitResultMap.end())
                continue;

            const bool isCurrentUpstream(LArClusterHelper::SortByPosition(pCurrentCluster, pTestCluster));

            // ATTN: Ensure that clusters are not contained within one another
            const float currentMinLayerZ(currentMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()), currentMaxLayerZ(currentMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ());
            const float testMinLayerZ(testMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()), testMaxLayerZ(testMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ());
            
            if (((currentMinLayerZ > testMinLayerZ) && (currentMaxLayerZ < testMaxLayerZ)) || ((testMinLayerZ > currentMinLayerZ) && (testMaxLayerZ < currentMaxLayerZ)))
                continue;

            CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
            if (!this->GetClusterMergingCoordinates(currentMicroFitIter->second, currentMacroFitIter->second, testMacroFitIter->second, currentMergePoint, currentMergeDirection, isCurrentUpstream) ||
                !this->GetClusterMergingCoordinates(testMicroFitIter->second, testMacroFitIter->second, currentMacroFitIter->second, testMergePoint, testMergeDirection, !isCurrentUpstream))
            {
                continue;
            }
            
            if ((isCurrentUpstream && !this->AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection)) ||
                (!isCurrentUpstream && !this->AreClustersAssociated(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection)))
            {
                continue;
            }

            if (isCurrentUpstream)
            {
                clusterAssociation = ClusterAssociation(pCurrentCluster, pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection);
            }
            else
            {
                clusterAssociation = ClusterAssociation(pTestCluster, pCurrentCluster, testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection);
            }

            foundAssociation = true;
            maxLength = lengthSum;
        }
    }

    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, CartesianVector &currentMergePosition, CartesianVector &currentMergeDirection, const bool isUpstream) const
{
    CartesianVector currentAverageDirection(0.f, 0.f, 0.f);
    CartesianVector associatedAverageDirection(0.f, 0.f, 0.f);
    currentMacroFitResult.GetGlobalDirection(currentMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), currentAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);    
    
    const LayerFitResultMap &currentMicroLayerFitResultMap(currentMicroFitResult.GetLayerFitResultMap());
    const int startLayer(isUpstream ? currentMicroFitResult.GetMaxLayer() : currentMicroFitResult.GetMinLayer());
    const int endLayer(isUpstream ? currentMicroFitResult.GetMinLayer() : currentMicroFitResult.GetMaxLayer());
    const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
    const int step(isUpstream ? -1 : 1);

    // ATTN: Search for stable region for which the local layer gradient agrees well with associated cluster global gradient
    unsigned int gradientStabilityHitWindow(std::ceil(currentMicroFitResult.GetCluster()->GetNCaloHits() * m_stableRegionClusterFraction)); // ISOBEL: SHOULD THIS BE A LAYER SPAN??
    unsigned int goodPositionCount(0);      

    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(currentMicroLayerFitResultMap.find(i));

        if (microIter == currentMicroLayerFitResultMap.end())
            continue;
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        currentMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
    
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            // ATTN: Set merge point and direction as that at the first layer in the stable region
            if (goodPositionCount == 0)
            {
                // ATTN: Cluster direction vectors must point to one another
                currentMergeDirection = currentAverageDirection * (isUpstream ? 1.f : -1.f);
                currentMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), currentMergePosition);
            }
            
            ++goodPositionCount;
        }
        else
        {
            goodPositionCount = 0;
        }

        if (goodPositionCount > gradientStabilityHitWindow)
            break;
        
        // ATTN: Abort merging process have not found a stable region 
        if (i == endLayer)
            return false;                                         
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::AreClustersAssociated(const CartesianVector &upstreamPoint, const CartesianVector &upstreamDirection, const CartesianVector &downstreamPoint,
    const CartesianVector &downstreamDirection) const
{
    if (downstreamPoint.GetZ() < upstreamPoint.GetZ())
        return false;

    const float separationDistance(std::sqrt(upstreamPoint.GetDistanceSquared(downstreamPoint)));
    if (separationDistance < m_minSeparationDistance)
        return false;
    
    // ATTN: Cluster directions point towards one another
    if (upstreamDirection.GetCosOpeningAngle(downstreamDirection * (-1.f)) < m_minDirectionDeviationCosAngle)
        return false;
    
    const CartesianVector predictedDownstreamPoint(upstreamPoint + (upstreamDirection * separationDistance));
    const float predictedDownstreamOffset((predictedDownstreamPoint - downstreamPoint).GetMagnitude());
    
    if (predictedDownstreamOffset > m_maxPredictedMergePointOffset)
        return false;

    const CartesianVector predictedUpstreamPoint(downstreamPoint + (downstreamDirection * separationDistance));
    const float predictedUpstreamOffset((predictedUpstreamPoint - upstreamPoint).GetMagnitude());

    if (predictedUpstreamOffset > m_maxPredictedMergePointOffset)
        return false;
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList,
    CaloHitVector &extrapolatedCaloHitVector, ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    const float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));

    // ATTN: Consider hits from merging clusters as merging point may not be cluster end but do not put them in clusterToCaloHitListMap
    for (const Cluster *const pCluster : *pClusterList) 
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;
                
                const float distanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(hitPosition - upstreamPoint).GetMagnitude());
                if (distanceFromLine > m_distanceFromLine)
                    continue;
                
                extrapolatedCaloHitVector.push_back(pCaloHit);

                if ((pCluster != clusterAssociation.GetUpstreamCluster()) && (pCluster != clusterAssociation.GetDownstreamCluster()))
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, CaloHitVector &extrapolatedCaloHitVector) const
{
    const CartesianVector &upstreamMergePoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamMergePoint(clusterAssociation.GetDownstreamMergePoint());
    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    const CartesianVector trackStep(trackDirection * m_lineSegmentLength);

    // To handle final segment - merge track remainder with preceding segment and if track remainder was more than half trackstep split into two
    const float trackLength((downstreamMergePoint - upstreamMergePoint).GetMagnitude());
    const unsigned int fullSegments(floor(trackLength / m_lineSegmentLength));
    const float lengthOfTrackRemainder(trackLength - (fullSegments * m_lineSegmentLength));

    // Sort hits by projected distance from the upstreamMergePoint
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), SortByDistanceAlongLine(upstreamMergePoint, trackDirection));

    unsigned int segmentsWithoutHits(0);
    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());

    for (unsigned int i = 0; i < (fullSegments + 1); ++i)
    {
        // if have run out of hits
        if (caloHitIter == extrapolatedCaloHitVector.end())
        {
            ++segmentsWithoutHits;
            
            if (segmentsWithoutHits > m_maxTrackGaps)
                return false;
            
            continue;
        }

        CartesianVector lowerBoundary(upstreamMergePoint + (trackStep * i));
        CartesianVector upperBoundary(upstreamMergePoint + (trackStep * (i + 1)));

        if (i >= fullSegments - 1)
        {
            if (lengthOfTrackRemainder > m_lineSegmentLength * 0.5f)
            {
                lowerBoundary = lowerBoundary - (trackStep * 0.5f * (i - fullSegments + 1.f)) + (trackDirection * 0.5f * (i - fullSegments + 1.f) * lengthOfTrackRemainder);
                upperBoundary = upperBoundary - (trackStep * 0.5f * (i - fullSegments + 2.f)) + (trackDirection * 0.5f * (i - fullSegments + 2.f) * lengthOfTrackRemainder);
            }
            else
            {
                upperBoundary = downstreamMergePoint;
            }
        }

        unsigned int hitsInSegment(0);
        while (this->IsInLineSegment(lowerBoundary, upperBoundary, (*caloHitIter)->GetPositionVector()))
        {
            ++hitsInSegment;
            ++caloHitIter;

            if (caloHitIter == extrapolatedCaloHitVector.end())
                break;
        }

        segmentsWithoutHits = hitsInSegment ? 0 : segmentsWithoutHits + 1;

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

        // leave loop early if final merged segment is not split into two
        if (i == (fullSegments - 1) && !(lengthOfTrackRemainder > m_lineSegmentLength * 0.5))
            return true;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float gradient = (-1.0)*(upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ() + gradient*upperBoundary.GetX())/gradient);
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ() + gradient*lowerBoundary.GetX())/gradient);

    const CartesianVector upper(xPointOnUpperLine, 0.f, point.GetZ());
    const CartesianVector lower(xPointOnLowerLine, 0.f, point.GetZ());

    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::CreateMainTrack(const ClusterAssociation &clusterAssociation, const CaloHitVector &extrapolatedCaloHitVector, const ClusterToCaloHitListMap &clusterToCaloHitListMap, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterVector &clusterVector) const
{
    ClusterList remnantClusterList;

    // Remove off axis hits from main track cluster pair
    const Cluster *const pMainTrackCluster(this->RemoveOffAxisHitsFromTrack(clusterAssociation.GetUpstreamCluster(), clusterAssociation.GetUpstreamMergePoint(),
        true, extrapolatedCaloHitVector, remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));
    const Cluster *const pClusterToDelete(this->RemoveOffAxisHitsFromTrack(clusterAssociation.GetDownstreamCluster(), clusterAssociation.GetDownstreamMergePoint(),
        false, extrapolatedCaloHitVector, remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));

    // Merge main track cluster pair together
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pClusterToDelete));

    // Determine the shower clusters that contain hits which belong to the main track
    ClusterVector showerClustersToFragment;
    for (const auto &entry : clusterToCaloHitListMap)
        showerClustersToFragment.push_back(entry.first);

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    // Remove hits from inbetween shower
    if (!showerClustersToFragment.empty())
    {
        for (const Cluster *const pShowerCluster : showerClustersToFragment)
        {
            this->UpdateForClusterDeletion(pShowerCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
            
            const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
            this->AddHitsToMainTrack(pShowerCluster, pMainTrackCluster, caloHitsToMerge, clusterAssociation, remnantClusterList);
        }
    }

    // Process the remnant clusters
    for (const Cluster *const pRemnantCluster : remnantClusterList)
    {
        if (pRemnantCluster->GetNCaloHits() == 1)
        {
            this->AddToNearestCluster(pRemnantCluster, pMainTrackCluster);
            continue;
        }
        
        if (this->IsClusterRemnantDisconnected(pRemnantCluster))
            this->FragmentRemnantCluster(pRemnantCluster);
    }

    this->UpdateAfterMainTrackCreation(pMainTrackCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackInEMShowerAlgorithm::RemoveOffAxisHitsFromTrack(const Cluster *const pCluster, const CartesianVector &splitPosition, const bool isUpstream,
    const CaloHitVector &extrapolatedCaloHitVector, ClusterList &remnantClusterList, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap,
    ClusterVector &clusterVector) const
{
    float rL(0.f), rT(0.f);
    const TwoDSlidingFitResult &microFitResult(microSlidingFitResultMap.at(pCluster));
    microFitResult.GetLocalPosition(splitPosition, rL, rT);

    const TwoDSlidingFitResult macroFitResult(macroSlidingFitResultMap.at(pCluster));
    CartesianVector averageDirection(0.f, 0.f, 0.f);
    macroFitResult.GetGlobalDirection(macroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), averageDirection);

    const float clusterGradient(averageDirection.GetZ()/averageDirection.GetX());           
    const float clusterIntercept(splitPosition.GetZ() - (clusterGradient * splitPosition.GetX()));

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pMainTrackCluster(nullptr), *pAboveCluster(nullptr), *pBelowCluster(nullptr);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            float thisL(0.f), thisT(0.f);

            microFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            const bool isAnExtrapolatedHit(std::find(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), pCaloHit) != extrapolatedCaloHitVector.end());
            const bool isAbove(hitPosition.GetZ() > ((clusterGradient * hitPosition.GetX()) + clusterIntercept));
            const bool isToRemove(!isAnExtrapolatedHit && (((thisL > rL) && isUpstream) || ((thisL < rL) && !isUpstream)));
                
            const Cluster *&pClusterToModify(isToRemove ? (isAbove ? pAboveCluster : pBelowCluster) : pMainTrackCluster);
            
            if (pClusterToModify)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));

                if (pClusterToModify != pMainTrackCluster)
                    remnantClusterList.push_back(pClusterToModify);  
            }
        }
    }

    this->UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
    
    // End fragmentation
    if (pAboveCluster || pBelowCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        return pMainTrackCluster;
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        return pCluster;
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddHitsToMainTrack(const Cluster *const pShowerCluster, const Cluster *const pMainTrackCluster, const CaloHitList &caloHitsToMerge,
    const ClusterAssociation &clusterAssociation, ClusterList &remnantClusterList) const
{
    if (pShowerCluster->GetNCaloHits() == caloHitsToMerge.size())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pShowerCluster));
        return;
    }

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pShowerCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pAboveCluster(nullptr), *pBelowCluster(nullptr);
    
    const float connectingLineGradient(clusterAssociation.GetConnectingLineDirection().GetZ()/clusterAssociation.GetConnectingLineDirection().GetX());         
    const float connectingLineIntercept(clusterAssociation.GetUpstreamMergePoint().GetZ() - (connectingLineGradient * clusterAssociation.GetUpstreamMergePoint().GetX()));

    const OrderedCaloHitList orderedCaloHitList(pShowerCluster->GetOrderedCaloHitList()); 
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const bool isAnExtrapolatedHit(std::find(caloHitsToMerge.begin(), caloHitsToMerge.end(), pCaloHit) != caloHitsToMerge.end());
            if (isAnExtrapolatedHit)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pShowerCluster, pCaloHit));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pMainTrackCluster, pCaloHit));
            }
            else
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                const bool isAbove(hitPosition.GetZ() > ((connectingLineGradient * hitPosition.GetX()) + connectingLineIntercept));
                const Cluster *&pClusterToModify(isAbove ? pAboveCluster : pBelowCluster);
            
                if (pClusterToModify)
                {
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
                }
                else
                {
                    PandoraContentApi::Cluster::Parameters parameters;
                    parameters.m_caloHitList.push_back(pCaloHit);
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));

                    remnantClusterList.push_back(pClusterToModify);
                }
            }
        }
    }

    // End fragmentation
    if (pAboveCluster || pBelowCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
    }
            
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsClusterRemnantDisconnected(const Cluster *const pRemnantCluster) const
{
    if (pRemnantCluster->GetNCaloHits() == 1)
        return false;

    unsigned int count(0);
    float maxSeparationDistance(-std::numeric_limits<float>::max());

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    CartesianVector previousHitPosition((*orderedCaloHitList.begin()->second->begin())->GetPositionVector());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            ++count;
            
            if (count == 1)
                continue;
            
            const float separationDistance(std::sqrt(pCaloHit->GetPositionVector().GetDistanceSquared(previousHitPosition)));

            if (separationDistance > maxSeparationDistance)
                maxSeparationDistance = separationDistance;

            previousHitPosition = pCaloHit->GetPositionVector();
        }
    }

    if (maxSeparationDistance > m_maxHitSeparationForConnectedCluster)
        return true;
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::FragmentRemnantCluster(const Cluster *const pRemnantCluster) const
{
    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pRemnantCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    ClusterList createdClusters;

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const Cluster *pClosestCluster(nullptr);
            
            if (!createdClusters.empty())
            {
                float closestDistance(std::numeric_limits<float>::max());

                for (const Cluster *const pCreatedCluster : createdClusters)
                {
                    const float separationDistance(LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pCreatedCluster));
                    if ((separationDistance < closestDistance) && (separationDistance < m_maxHitDistanceFromCluster))
                    {
                        closestDistance = separationDistance;
                        pClosestCluster = pCreatedCluster;
                    }
                }
            }

            if (pClosestCluster)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClosestCluster, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClosestCluster));
                createdClusters.push_back(pClosestCluster);
            }
        }
    }

    // End fragmentation
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName)); 
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddToNearestCluster(const Cluster *const pClusterToMerge, const Cluster *const pClusterToEnlarge) const
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    float closestDistance(std::numeric_limits<float>::max());
    const Cluster *pClosestCluster(nullptr);

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pClusterToMerge)
            continue;
        
        float separationDistance(LArClusterHelper::GetClosestDistance(pClusterToMerge, pCluster));

        if (separationDistance < closestDistance)
        {
            if ((pCluster == pClusterToEnlarge) && (separationDistance > m_maxDistanceFromMainTrack))
                continue;
            
            pClosestCluster = pCluster;
            closestDistance = separationDistance;
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClosestCluster, pClusterToMerge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateForClusterDeletion(const Cluster *const pCluster, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    std::vector<TwoDSlidingFitResultMap*> slidingFitResultMapVector({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    
    for (TwoDSlidingFitResultMap *const pSlidingFitResultMap : slidingFitResultMapVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitToDelete(pSlidingFitResultMap->find(pCluster));
        if (fitToDelete != pSlidingFitResultMap->end())
            pSlidingFitResultMap->erase(fitToDelete);
    }

    ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pCluster));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateAfterMainTrackCreation(const Cluster *const pMainTrackCluster, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    clusterVector.push_back(pMainTrackCluster);
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    try
    {
        TwoDSlidingFitResult microSlidingFitResult(pMainTrackCluster, m_slidingFitWindow, slidingFitPitch);
        TwoDSlidingFitResult macroSlidingFitResult(pMainTrackCluster, m_macroSlidingFitWindow, slidingFitPitch);

        (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pMainTrackCluster, microSlidingFitResult));
        (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pMainTrackCluster, macroSlidingFitResult));
    }
    catch (StatusCodeException &) {}
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackInEMShowerAlgorithm::ClusterAssociation::ClusterAssociation() :
    m_pUpstreamCluster(nullptr),
    m_pDownstreamCluster(nullptr),
    m_upstreamMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_upstreamMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_connectingLineDirection(CartesianVector(0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackInEMShowerAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pUpstreamCluster, const Cluster *const pDownstreamCluster,const CartesianVector &upstreamMergePoint,
        const CartesianVector &upstreamMergeDirection, const CartesianVector &downstreamMergePoint, const CartesianVector &downstreamMergeDirection) :
    m_pUpstreamCluster(pUpstreamCluster),
    m_pDownstreamCluster(pDownstreamCluster),
    m_upstreamMergePoint(upstreamMergePoint),
    m_upstreamMergeDirection(upstreamMergeDirection),
    m_downstreamMergePoint(downstreamMergePoint),
    m_downstreamMergeDirection(downstreamMergeDirection),
    m_connectingLineDirection(0.f, 0.f, 0.f)
{
    const CartesianVector connectingLineDirection(m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::SortByDistanceAlongLine::operator() (const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs)
{
    const CartesianVector lhsDistanceVector(pLhs->GetPositionVector() - m_startPoint);
    const CartesianVector rhsDistanceVector(pRhs->GetPositionVector() - m_startPoint);

    const float lhsProjectedDistance(lhsDistanceVector.GetDotProduct(m_lineDirection));
    const float rhsProjectedDistance(rhsDistanceVector.GetDotProduct(m_lineDirection));

    return (lhsProjectedDistance < rhsProjectedDistance);
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode TrackInEMShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPredictedMergePointOffset", m_maxPredictedMergePointOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparationDistance", m_minSeparationDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergePointMinCosAngleDeviation", m_mergePointMinCosAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLengthSum", m_minClusterLengthSum));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFromLine", m_distanceFromLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackGaps", m_maxTrackGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LineSegmentLength", m_lineSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitSeparationForConnectedCluster", m_maxHitSeparationForConnectedCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromMainTrack", m_maxDistanceFromMainTrack));        

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMainLoopIterations", m_maxMainLoopIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MacroSlidingFitWindow", m_macroSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StableRegionClusterFraction", m_stableRegionClusterFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
