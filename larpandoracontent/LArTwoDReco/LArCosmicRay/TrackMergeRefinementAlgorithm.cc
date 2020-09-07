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
    m_maxLoopIterations(10),    
    m_minClusterLengthSum(75.f),
    m_minSeparationDistance(0.f),
    m_minDirectionDeviationCosAngle(0.99f),
    m_maxPredictedMergePointOffset(5.f),
    m_distanceToLine(0.35f),    
    m_boundaryTolerance(2.f)  
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
StatusCode TrackMergeRefinementAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    ClusterVector clusterVector;
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    SlidingFitResultMapPair slidingFitResultMapPair({&microSlidingFitResultMap, &macroSlidingFitResultMap});

    this->InitialiseContainers(pClusterList, LArClusterHelper::SortByNHits, clusterVector, slidingFitResultMapPair);

    // ATTN: Keep track of created main track clusters so their hits can be protected in future iterations    
    unsigned int loopIterations(0);
    ClusterList createdMainTrackClusters;    
    while(loopIterations < m_maxLoopIterations) 
    {
        ++loopIterations;

        ClusterPairAssociation clusterAssociation;
        if(!this->FindBestClusterAssociation(clusterVector, slidingFitResultMapPair, clusterAssociation))
            break;

        ClusterList unavailableProtectedClusters;
        this->GetUnavailableProtectedClusters(clusterAssociation, createdMainTrackClusters, unavailableProtectedClusters);
        
        ClusterToCaloHitListMap clusterToCaloHitListMap;
        this->GetHitsInBoundingBox(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetDownstreamMergePoint(), pClusterList, clusterToCaloHitListMap,
            unavailableProtectedClusters, m_distanceToLine);

        if (!this->AreExtrapolatedHitsGood(clusterToCaloHitListMap, clusterAssociation))
        {
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
    
bool TrackMergeRefinementAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    ClusterPairAssociation &clusterAssociation) const
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

            if ((isCurrentUpstream && !this->AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection)) ||
                (!isCurrentUpstream && !this->AreClustersAssociated(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection)))
            {
                continue;
            }


            if (isCurrentUpstream)
            {
                //this->RepositionIfInGap(currentMergeDirection, testMergeDirection, currentMergePoint, testMergePoint);
                clusterAssociation = ClusterPairAssociation(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection * (-1.f), pCurrentCluster, pTestCluster);
            }
            else
            {
                //this->RepositionIfInGap(testMergeDirection, currentMergeDirection, testMergePoint, currentMergePoint);
                clusterAssociation = ClusterPairAssociation(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection * (-1.f), pTestCluster, pCurrentCluster);
            }

            foundAssociation = true;
            maxLength = lengthSum;
        }
    }
    
    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackMergeRefinementAlgorithm::RepositionIfInGap(const CartesianVector &upstreamMergeDirection, const CartesianVector &downstreamMergeDirection,
    CartesianVector &upstreamMergePoint, CartesianVector &downstreamMergePoint) const
{
    CartesianVector &lowerXMergePoint(upstreamMergePoint.GetX() < downstreamMergePoint.GetX() ? upstreamMergePoint : downstreamMergePoint);
    const float lowerGradient(upstreamMergePoint.GetX() < downstreamMergePoint.GetX() ?
        upstreamMergeDirection.GetZ() / upstreamMergeDirection.GetX() : downstreamMergeDirection.GetZ() / downstreamMergeDirection.GetX());                                      
    CartesianVector &higherXMergePoint(upstreamMergePoint.GetX() < downstreamMergePoint.GetX() ? downstreamMergePoint : upstreamMergePoint);
    const float higherGradient(upstreamMergePoint.GetX() < downstreamMergePoint.GetX() ?
        downstreamMergeDirection.GetZ() / downstreamMergeDirection.GetX() : upstreamMergeDirection.GetZ() / upstreamMergeDirection.GetX());
    const float upstreamGradient(upstreamMergeDirection.GetZ() / upstreamMergeDirection.GetX());
    const float downstreamGradient(downstreamMergeDirection.GetZ() / downstreamMergeDirection.GetX());    

    DetectorGapList detectorGapList(this->GetPandora().GetGeometry()->GetDetectorGapList());
    for (const DetectorGap *const pDetectorGap : detectorGapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pDetectorGap));

        if (pLineGap)
        {
            LineGapType lineGapType(pLineGap->GetLineGapType());
            
            if (lineGapType == TPC_DRIFT_GAP)
            {
                if ((pLineGap->GetLineStartX() < lowerXMergePoint.GetX()) && (pLineGap->GetLineEndX() > lowerXMergePoint.GetX()))
                    lowerXMergePoint = CartesianVector(pLineGap->GetLineEndX(), 0.f, lowerXMergePoint.GetZ() + lowerGradient * (pLineGap->GetLineEndX() - lowerXMergePoint.GetX()));

                if ((pLineGap->GetLineStartX() < higherXMergePoint.GetX()) && (pLineGap->GetLineEndX() > higherXMergePoint.GetX()))
                    higherXMergePoint = CartesianVector(pLineGap->GetLineStartX(), 0.f, higherXMergePoint.GetZ() + higherGradient * (pLineGap->GetLineStartX() - higherXMergePoint.GetX()));
            }

            if ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W))
            {
                if ((pLineGap->GetLineStartZ() < upstreamMergePoint.GetZ()) && (pLineGap->GetLineEndZ() > upstreamMergePoint.GetZ()))
                    upstreamMergePoint = CartesianVector((pLineGap->GetLineEndZ() - upstreamMergePoint.GetZ()) / upstreamGradient, 0.f, pLineGap->GetLineEndZ());

                if ((pLineGap->GetLineStartZ() < downstreamMergePoint.GetZ()) && (pLineGap->GetLineEndZ() > downstreamMergePoint.GetZ()))
                    downstreamMergePoint = CartesianVector((pLineGap->GetLineStartZ() - downstreamMergePoint.GetZ()) / downstreamGradient, 0.f, pLineGap->GetLineStartZ());
            }
        }
    }
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

void TrackMergeRefinementAlgorithm::GetUnavailableProtectedClusters(const ClusterPairAssociation &clusterAssociation, const ClusterList &createdMainTrackClusters,
    ClusterList &protectedClusters) const
{
    for (const Cluster *const pMainTrackCluster : createdMainTrackClusters)
    {
        if ((pMainTrackCluster != clusterAssociation.GetUpstreamCluster()) && (pMainTrackCluster != clusterAssociation.GetDownstreamCluster()))
            protectedClusters.push_back(pMainTrackCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackMergeRefinementAlgorithm::AreExtrapolatedHitsNearBoundaries(const CaloHitVector &extrapolatedHitVector, ClusterAssociation &clusterAssociation) const
{
    if (extrapolatedHitVector.empty())
    {
        const float separationDistance((clusterAssociation.GetUpstreamMergePoint() - clusterAssociation.GetDownstreamMergePoint()).GetMagnitude());
        return (separationDistance < m_boundaryTolerance);
    }

    if (!this->IsNearBoundary(extrapolatedHitVector.front(), clusterAssociation.GetUpstreamMergePoint(), m_boundaryTolerance))
        return false;

    if (!this->IsNearBoundary(extrapolatedHitVector.back(), clusterAssociation.GetDownstreamMergePoint(), m_boundaryTolerance))
        return false;
    
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

    // ATTN: Cleanup containers - choose to add created clusters back into containers
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
        "MaxLoopIterations", m_maxLoopIterations));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLengthSum", m_minClusterLengthSum));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparationDistance", m_minSeparationDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPredictedMergePointOffset", m_maxPredictedMergePointOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceToLine", m_distanceToLine));       

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BoundaryTolerance", m_boundaryTolerance));

    return TrackRefinementBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------     

}// namespace lar_content
