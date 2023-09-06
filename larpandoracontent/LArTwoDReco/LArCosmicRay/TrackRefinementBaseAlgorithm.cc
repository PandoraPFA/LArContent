/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackRefinementBaseAlgorithm.cc
 *
 *  @brief  Implementation of the track refinement base class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackRefinementBaseAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"

using namespace pandora;

namespace lar_content
{

TrackRefinementBaseAlgorithm::TrackRefinementBaseAlgorithm() :
    m_minClusterLength(15.f),
    m_microSlidingFitWindow(20),
    m_macroSlidingFitWindow(1000),
    m_stableRegionClusterFraction(0.05f),
    m_mergePointMinCosAngleDeviation(0.999f),
    m_minHitFractionForHitRemoval(0.05f),
    m_maxDistanceFromMainTrack(0.75f),
    m_maxHitDistanceFromCluster(4.f),
    m_maxHitSeparationForConnectedCluster(4.f),
    m_maxTrackGaps(3),
    m_lineSegmentLength(3.f),
    m_hitWidthMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TrackRefinementBaseAlgorithm::InitialiseContainers(const ClusterList *pClusterList, const T sortFunction, ClusterVector &clusterVector,
    SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
            const TwoDSlidingFitResult microSlidingFitResult(pCluster, m_microSlidingFitWindow, slidingFitPitch);
            const TwoDSlidingFitResult macroSlidingFitResult(pCluster, m_macroSlidingFitWindow, slidingFitPitch);

            slidingFitResultMapPair.first->insert(TwoDSlidingFitResultMap::value_type(pCluster, microSlidingFitResult));
            slidingFitResultMapPair.second->insert(TwoDSlidingFitResultMap::value_type(pCluster, macroSlidingFitResult));
            clusterVector.push_back(pCluster);
        }
        catch (const StatusCodeException &)
        {
        }
    }

    std::sort(clusterVector.begin(), clusterVector.end(), sortFunction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &clusterMicroFitResult,
    const TwoDSlidingFitResult &clusterMacroFitResult, const TwoDSlidingFitResult &associatedMacroFitResult, const bool isEndUpstream,
    CartesianVector &clusterMergePosition, CartesianVector &clusterMergeDirection) const
{
    CartesianVector clusterAverageDirection(0.f, 0.f, 0.f), associatedAverageDirection(0.f, 0.f, 0.f);
    clusterMacroFitResult.GetGlobalDirection(clusterMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);

    const LayerFitResultMap &clusterMicroLayerFitResultMap(clusterMicroFitResult.GetLayerFitResultMap());
    const int startLayer(isEndUpstream ? clusterMicroFitResult.GetMinLayer() : clusterMicroFitResult.GetMaxLayer());
    const int endLayer(isEndUpstream ? clusterMicroFitResult.GetMaxLayer() : clusterMicroFitResult.GetMinLayer());
    const int loopTerminationLayer(endLayer + (isEndUpstream ? 1 : -1));
    const int step(isEndUpstream ? 1 : -1);

    // ATTN: Search for stable region for which the local layer gradient agrees well with associated cluster global gradient
    unsigned int gradientStabilityWindow(std::ceil(clusterMicroLayerFitResultMap.size() * m_stableRegionClusterFraction));
    unsigned int goodLayerCount(0);

    clusterMicroLayerFitResultMap.at(endLayer);

    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        clusterMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);

        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            // ATTN: Set merge point and direction as that at the first layer in the stable region
            if (goodLayerCount == 0)
            {
                clusterMergeDirection = clusterAverageDirection;
                PANDORA_THROW_RESULT_IF(
                    STATUS_CODE_SUCCESS, !=, clusterMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), clusterMergePosition));
            }

            ++goodLayerCount;
        }
        else
        {
            goodLayerCount = 0;
        }

        if (goodLayerCount > gradientStabilityWindow)
            break;

        // ATTN: Abort merging process have not found a stable region
        if (i == endLayer)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRefinementBaseAlgorithm::GetHitsInBoundingBox(const CartesianVector &firstCorner, const CartesianVector &secondCorner,
    const ClusterList *const pClusterList, ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList &unavailableProtectedClusters, const float distanceToLine) const
{
    const float minX(std::min(firstCorner.GetX(), secondCorner.GetX())), maxX(std::max(firstCorner.GetX(), secondCorner.GetX()));
    const float minZ(std::min(firstCorner.GetZ(), secondCorner.GetZ())), maxZ(std::max(firstCorner.GetZ(), secondCorner.GetZ()));

    CartesianVector connectingLineDirection(firstCorner - secondCorner);
    connectingLineDirection = connectingLineDirection.GetUnitVector();

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (std::find(unavailableProtectedClusters.begin(), unavailableProtectedClusters.end(), pCluster) != unavailableProtectedClusters.end())
            continue;

        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                CartesianVector hitPosition(m_hitWidthMode ? LArHitWidthHelper::GetClosestPointToLine2D(firstCorner, connectingLineDirection, pCaloHit)
                                                           : pCaloHit->GetPositionVector());

                if (!this->IsInBoundingBox(minX, maxX, minZ, maxZ, hitPosition))
                    continue;

                if (distanceToLine > 0.f)
                {
                    if (!this->IsCloseToLine(hitPosition, firstCorner, connectingLineDirection, distanceToLine))
                        continue;
                }

                clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::IsInBoundingBox(
    const float minX, const float maxX, const float minZ, const float maxZ, const CartesianVector &hitPosition) const
{
    return !((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < minZ) || (hitPosition.GetZ() > maxZ));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart,
    const CartesianVector &lineDirection, const float distanceToLine) const
{
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());

    return (transverseDistanceFromLine < distanceToLine);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::AreExtrapolatedHitsGood(const ClusterToCaloHitListMap &clusterToCaloHitListMap, ClusterAssociation &clusterAssociation) const
{
    CaloHitVector extrapolatedHitVector;
    for (const auto &entry : clusterToCaloHitListMap)
        extrapolatedHitVector.insert(extrapolatedHitVector.begin(), entry.second.begin(), entry.second.end());

    // ATTN: Extrapolated hit checks require extrapolatedHitVector to be ordered from upstream -> downstream merge point
    std::sort(extrapolatedHitVector.begin(), extrapolatedHitVector.end(),
        SortByDistanceAlongLine(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetConnectingLineDirection(), m_hitWidthMode));

    if (!this->AreExtrapolatedHitsNearBoundaries(extrapolatedHitVector, clusterAssociation))
        return false;

    if (clusterToCaloHitListMap.empty())
        return true;

    if (!this->IsTrackContinuous(clusterAssociation, extrapolatedHitVector))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::IsNearBoundary(const CaloHit *const pCaloHit, const CartesianVector &boundaryPosition2D, const float boundaryTolerance) const
{
    const float distanceToBoundary(LArHitWidthHelper::GetClosestDistanceToPoint2D(pCaloHit, boundaryPosition2D));

    return (distanceToBoundary < boundaryTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, const CaloHitVector &extrapolatedCaloHitVector) const
{
    CartesianPointVector trackSegmentBoundaries;
    this->GetTrackSegmentBoundaries(clusterAssociation, trackSegmentBoundaries);

    if (trackSegmentBoundaries.size() < 2)
    {
        std::cout << "TrackInEMShowerAlgorithm: Less than two track segment boundaries" << std::endl;
        throw STATUS_CODE_NOT_ALLOWED;
    }

    unsigned int segmentsWithoutHits(0);
    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());
    CartesianVector hitPosition(m_hitWidthMode ? LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(),
                                                     clusterAssociation.GetConnectingLineDirection(), *caloHitIter)
                                               : (*caloHitIter)->GetPositionVector());

    for (unsigned int i = 0; i < (trackSegmentBoundaries.size() - 1); ++i)
    {
        if (caloHitIter == extrapolatedCaloHitVector.end())
        {
            ++segmentsWithoutHits;

            if (segmentsWithoutHits > m_maxTrackGaps)
                return false;

            continue;
        }

        unsigned int hitsInSegment(0);
        while (this->IsInLineSegment(trackSegmentBoundaries.at(i), trackSegmentBoundaries.at(i + 1), hitPosition))
        {
            ++hitsInSegment;
            ++caloHitIter;

            if (caloHitIter == extrapolatedCaloHitVector.end())
                break;

            hitPosition = m_hitWidthMode ? LArHitWidthHelper::GetClosestPointToLine2D(clusterAssociation.GetUpstreamMergePoint(),
                                               clusterAssociation.GetConnectingLineDirection(), *caloHitIter)
                                         : (*caloHitIter)->GetPositionVector();
        }

        segmentsWithoutHits = hitsInSegment ? 0 : segmentsWithoutHits + 1;

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRefinementBaseAlgorithm::GetTrackSegmentBoundaries(const ClusterAssociation &clusterAssociation, CartesianPointVector &trackSegmentBoundaries) const
{
    if (m_lineSegmentLength < std::numeric_limits<float>::epsilon())
    {
        std::cout << "TrackInEMShowerAlgorithm: Line segment length must be positive and nonzero" << std::endl;
        throw STATUS_CODE_INVALID_PARAMETER;
    }

    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    float trackLength((clusterAssociation.GetUpstreamMergePoint() - clusterAssociation.GetDownstreamMergePoint()).GetMagnitude());

    // ATTN: Consider the existence of gaps where no hits will be found
    DetectorGapList consideredGaps;
    trackLength -= this->DistanceInGap(
        clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetDownstreamMergePoint(), trackDirection, consideredGaps);
    consideredGaps.clear();

    const int fullSegments(std::floor(trackLength / m_lineSegmentLength));

    if (fullSegments == 0)
    {
        trackSegmentBoundaries = {clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetDownstreamMergePoint()};
        return;
    }

    // ATTN: To handle final segment merge track remainder with preceding segment and if track remainder was more than half of the segment length split into two
    const float lengthOfTrackRemainder(trackLength - (fullSegments * m_lineSegmentLength));
    const bool splitFinalSegment(lengthOfTrackRemainder > m_lineSegmentLength * 0.5f);
    const int numberOfBoundaries(fullSegments + (splitFinalSegment ? 2 : 1));

    CartesianVector segmentUpstreamBoundary(clusterAssociation.GetUpstreamMergePoint()), segmentDownstreamBoundary(segmentUpstreamBoundary);
    this->RepositionIfInGap(trackDirection, segmentUpstreamBoundary);
    trackSegmentBoundaries.push_back(segmentUpstreamBoundary);

    for (int i = 1; i < numberOfBoundaries; ++i)
    {
        if (i < fullSegments)
        {
            segmentDownstreamBoundary += trackDirection * m_lineSegmentLength;
        }
        else
        {
            if (splitFinalSegment)
            {
                segmentDownstreamBoundary += trackDirection * (m_lineSegmentLength + lengthOfTrackRemainder) * 0.5f;
            }
            else
            {
                segmentDownstreamBoundary += trackDirection * (m_lineSegmentLength + lengthOfTrackRemainder);
            }
        }

        float distanceInGap(this->DistanceInGap(segmentUpstreamBoundary, segmentDownstreamBoundary, trackDirection, consideredGaps));
        while (distanceInGap > std::numeric_limits<float>::epsilon())
        {
            this->RepositionIfInGap(trackDirection, segmentDownstreamBoundary);
            segmentDownstreamBoundary += trackDirection * distanceInGap;
            distanceInGap = this->DistanceInGap(segmentUpstreamBoundary, segmentDownstreamBoundary, trackDirection, consideredGaps);
        }

        trackSegmentBoundaries.push_back(segmentDownstreamBoundary);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRefinementBaseAlgorithm::RepositionIfInGap(const CartesianVector &mergeDirection, CartesianVector &trackPoint) const
{
    const DetectorGapList detectorGapList(this->GetPandora().GetGeometry()->GetDetectorGapList());
    for (const DetectorGap *const pDetectorGap : detectorGapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pDetectorGap));

        if (pLineGap)
        {
            const LineGapType lineGapType(pLineGap->GetLineGapType());

            if (lineGapType == TPC_DRIFT_GAP)
            {
                if ((pLineGap->GetLineStartX() < trackPoint.GetX()) && (pLineGap->GetLineEndX() > trackPoint.GetX()))
                {
                    if (std::fabs(mergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
                        throw StatusCodeException(STATUS_CODE_FAILURE);

                    const float gradient(mergeDirection.GetZ() / mergeDirection.GetX());

                    trackPoint = CartesianVector(
                        pLineGap->GetLineEndX(), 0.f, trackPoint.GetZ() + (gradient * (pLineGap->GetLineEndX() - trackPoint.GetX())));
                }
            }

            if ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W))
            {
                if ((pLineGap->GetLineStartZ() < trackPoint.GetZ()) && (pLineGap->GetLineEndZ() > trackPoint.GetZ()))
                {
                    if (std::fabs(mergeDirection.GetZ()) < std::numeric_limits<float>::epsilon())
                        throw StatusCodeException(STATUS_CODE_FAILURE);

                    const float gradient(mergeDirection.GetX() / mergeDirection.GetZ());

                    trackPoint = CartesianVector(
                        trackPoint.GetX() + (gradient * (pLineGap->GetLineEndZ() - trackPoint.GetZ())), 0.f, pLineGap->GetLineEndZ());
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackRefinementBaseAlgorithm::DistanceInGap(const CartesianVector &upstreamPoint, const CartesianVector &downstreamPoint,
    const CartesianVector &connectingLine, DetectorGapList &consideredGaps) const
{
    const CartesianVector &lowerXPoint(upstreamPoint.GetX() < downstreamPoint.GetX() ? upstreamPoint : downstreamPoint);
    const CartesianVector &higherXPoint(upstreamPoint.GetX() < downstreamPoint.GetX() ? downstreamPoint : upstreamPoint);

    const float cosAngleToX(std::fabs(connectingLine.GetDotProduct(CartesianVector(1.f, 0.f, 0.f))));
    const float cosAngleToZ(std::fabs(connectingLine.GetDotProduct(CartesianVector(0.f, 0.f, 1.f))));

    float distanceInGaps(0.f);
    const DetectorGapList detectorGapList(this->GetPandora().GetGeometry()->GetDetectorGapList());
    for (const DetectorGap *const pDetectorGap : detectorGapList)
    {
        if (std::find(consideredGaps.begin(), consideredGaps.end(), pDetectorGap) != consideredGaps.end())
            continue;

        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pDetectorGap));

        if (pLineGap)
        {
            const LineGapType lineGapType(pLineGap->GetLineGapType());

            if (lineGapType == TPC_DRIFT_GAP)
            {
                float xDistanceInGap(0.f);

                if ((pLineGap->GetLineStartX() > lowerXPoint.GetX()) && (pLineGap->GetLineEndX() < higherXPoint.GetX()))
                {
                    xDistanceInGap = (pLineGap->GetLineEndX() - pLineGap->GetLineStartX());
                }
                else if ((pLineGap->GetLineStartX() < lowerXPoint.GetX()) && (pLineGap->GetLineEndX() > lowerXPoint.GetX()))
                {
                    xDistanceInGap = (pLineGap->GetLineEndX() - lowerXPoint.GetX());
                }
                else if ((pLineGap->GetLineStartX() < higherXPoint.GetX()) && (pLineGap->GetLineEndX() > higherXPoint.GetX()))
                {
                    xDistanceInGap = (higherXPoint.GetX() - pLineGap->GetLineStartX());
                }
                else
                {
                    continue;
                }

                if (std::fabs(cosAngleToX) < std::numeric_limits<float>::epsilon())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                distanceInGaps += (xDistanceInGap / cosAngleToX);

                consideredGaps.push_back(pDetectorGap);
            }

            if ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W))
            {
                float zDistanceInGap(0.f);

                if ((pLineGap->GetLineStartZ() > upstreamPoint.GetZ()) && (pLineGap->GetLineEndZ() < downstreamPoint.GetZ()))
                {
                    zDistanceInGap = pLineGap->GetLineEndZ() - pLineGap->GetLineStartZ();
                }
                else if ((pLineGap->GetLineStartZ() < upstreamPoint.GetZ()) && (pLineGap->GetLineEndZ() > upstreamPoint.GetZ()))
                {
                    zDistanceInGap = pLineGap->GetLineEndZ() - upstreamPoint.GetZ();
                }
                else if ((pLineGap->GetLineStartZ() < downstreamPoint.GetZ()) && (pLineGap->GetLineEndZ() > downstreamPoint.GetZ()))
                {
                    zDistanceInGap = downstreamPoint.GetZ() - pLineGap->GetLineStartZ();
                }
                else
                {
                    continue;
                }

                if (std::fabs(cosAngleToZ) < std::numeric_limits<float>::epsilon())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                distanceInGaps += (zDistanceInGap / cosAngleToZ);

                consideredGaps.push_back(pDetectorGap);
            }
        }
    }

    return distanceInGaps;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::IsInLineSegment(
    const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float segmentBoundaryGradient = (-1.f) * (upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ()) / segmentBoundaryGradient + upperBoundary.GetX());
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ()) / segmentBoundaryGradient + lowerBoundary.GetX());

    if (std::fabs(xPointOnUpperLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;

    if (std::fabs(xPointOnLowerLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;

    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackRefinementBaseAlgorithm::RemoveOffAxisHitsFromTrack(const Cluster *const pCluster, const CartesianVector &splitPosition,
    const bool isEndUpstream, const ClusterToCaloHitListMap &clusterToCaloHitListMap, ClusterList &remnantClusterList,
    TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    float rL(0.f), rT(0.f);
    const TwoDSlidingFitResult &microFitResult(microSlidingFitResultMap.at(pCluster));
    microFitResult.GetLocalPosition(splitPosition, rL, rT);

    const TwoDSlidingFitResult &macroFitResult(macroSlidingFitResultMap.at(pCluster));
    CartesianVector averageDirection(0.f, 0.f, 0.f);
    macroFitResult.GetGlobalDirection(macroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), averageDirection);

    const bool isVertical(std::fabs(averageDirection.GetX()) < std::numeric_limits<float>::epsilon());
    const float clusterGradient(isVertical ? 0.f : averageDirection.GetZ() / averageDirection.GetX());
    const float clusterIntercept(isVertical ? splitPosition.GetX() : splitPosition.GetZ() - (clusterGradient * splitPosition.GetX()));

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pCluster);
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pMainTrackCluster(nullptr), *pAboveCluster(nullptr), *pBelowCluster(nullptr);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            float thisL(0.f), thisT(0.f);

            microFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            bool isAnExtrapolatedHit(false);
            const auto extrapolatedCaloHitIter(clusterToCaloHitListMap.find(pCluster));

            if (extrapolatedCaloHitIter != clusterToCaloHitListMap.end())
                isAnExtrapolatedHit = std::find(extrapolatedCaloHitIter->second.begin(), extrapolatedCaloHitIter->second.end(), pCaloHit) !=
                                      extrapolatedCaloHitIter->second.end();

            const bool isAbove(((clusterGradient * hitPosition.GetX()) + clusterIntercept) < (isVertical ? hitPosition.GetX() : hitPosition.GetZ()));
            const bool isToRemove(!isAnExtrapolatedHit && (((thisL < rL) && isEndUpstream) || ((thisL > rL) && !isEndUpstream)));

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

void TrackRefinementBaseAlgorithm::AddHitsToMainTrack(const Cluster *const pMainTrackCluster, const Cluster *const pShowerCluster,
    const CaloHitList &caloHitsToMerge, const ClusterAssociation &clusterAssociation, ClusterList &remnantClusterList) const
{
    // To ignore crossing CR muon or test beam tracks
    if (((static_cast<float>(caloHitsToMerge.size()) / static_cast<float>(pShowerCluster->GetNCaloHits())) < m_minHitFractionForHitRemoval) &&
        (LArClusterHelper::GetLengthSquared(pShowerCluster) > m_minClusterLength))
        return;

    if (pShowerCluster->GetNCaloHits() == caloHitsToMerge.size())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pShowerCluster));
        return;
    }

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pShowerCluster);
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pAboveCluster(nullptr), *pBelowCluster(nullptr);

    const bool isVertical(std::fabs(clusterAssociation.GetConnectingLineDirection().GetX()) < std::numeric_limits<float>::epsilon());
    const float connectingLineGradient(
        isVertical ? 0.f : clusterAssociation.GetConnectingLineDirection().GetZ() / clusterAssociation.GetConnectingLineDirection().GetX());
    const float connectingLineIntercept(isVertical ? clusterAssociation.GetUpstreamMergePoint().GetX()
                                                   : clusterAssociation.GetUpstreamMergePoint().GetZ() -
                                                         (connectingLineGradient * clusterAssociation.GetUpstreamMergePoint().GetX()));

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
                const bool isAbove(((connectingLineGradient * hitPosition.GetX()) + connectingLineIntercept) <
                                   (isVertical ? hitPosition.GetX() : hitPosition.GetZ()));
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

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRefinementBaseAlgorithm::ProcessRemnantClusters(const ClusterList &remnantClusterList, const Cluster *const pMainTrackCluster,
    const ClusterList *const pClusterList, ClusterList &createdClusters) const
{
    ClusterList fragmentedClusters;
    for (const Cluster *const pRemnantCluster : remnantClusterList)
    {
        if (this->IsClusterRemnantDisconnected(pRemnantCluster))
        {
            this->FragmentRemnantCluster(pRemnantCluster, fragmentedClusters);
        }
        else
        {
            fragmentedClusters.push_back(pRemnantCluster);
        }
    }

    for (const Cluster *const pFragmentedCluster : fragmentedClusters)
    {
        if ((pFragmentedCluster->GetNCaloHits() == 1) && (this->AddToNearestCluster(pFragmentedCluster, pMainTrackCluster, pClusterList)))
            continue;

        createdClusters.push_back(pFragmentedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::AddToNearestCluster(
    const Cluster *const pClusterToMerge, const Cluster *const pMainTrackCluster, const ClusterList *const pClusterList) const
{
    const Cluster *pClosestCluster(nullptr);
    float closestDistance(std::numeric_limits<float>::max());

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pClusterToMerge)
            continue;

        const float separationDistance(LArClusterHelper::GetClosestDistance(pClusterToMerge, pCluster));

        if (separationDistance < closestDistance)
        {
            if ((pCluster == pMainTrackCluster) && (separationDistance > m_maxDistanceFromMainTrack))
                continue;

            pClosestCluster = pCluster;
            closestDistance = separationDistance;
        }
    }

    if (closestDistance < m_maxHitDistanceFromCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClosestCluster, pClusterToMerge));
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::IsClusterRemnantDisconnected(const Cluster *const pRemnantCluster) const
{
    if (pRemnantCluster->GetNCaloHits() == 1)
        return false;

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    const CaloHit *pPreviousCaloHit(orderedCaloHitList.begin()->second->front());

    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            if (pCaloHit == pPreviousCaloHit)
                continue;

            const float separationDistanceSquared(pCaloHit->GetPositionVector().GetDistanceSquared(pPreviousCaloHit->GetPositionVector()));

            if (separationDistanceSquared > (m_maxHitSeparationForConnectedCluster * m_maxHitSeparationForConnectedCluster))
                return true;

            pPreviousCaloHit = pCaloHit;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRefinementBaseAlgorithm::FragmentRemnantCluster(const Cluster *const pRemnantCluster, ClusterList &fragmentedClusterList) const
{
    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pRemnantCluster);
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

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
    if (createdClusters.empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        fragmentedClusterList.push_back(pRemnantCluster);
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        fragmentedClusterList.insert(fragmentedClusterList.begin(), createdClusters.begin(), createdClusters.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TrackRefinementBaseAlgorithm::UpdateContainers(const ClusterList &clustersToAdd, const ClusterList &clustersToDelete,
    const T sortFunction, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    //ATTN: Very important to first delete pointers from containers
    for (const Cluster *const pClusterToDelete : clustersToDelete)
        this->RemoveClusterFromContainers(pClusterToDelete, clusterVector, slidingFitResultMapPair);

    this->InitialiseContainers(&clustersToAdd, sortFunction, clusterVector, slidingFitResultMapPair);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRefinementBaseAlgorithm::RemoveClusterFromContainers(
    const Cluster *const pClusterToRemove, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    const TwoDSlidingFitResultMap::const_iterator microFitToDelete(slidingFitResultMapPair.first->find(pClusterToRemove));
    if (microFitToDelete != slidingFitResultMapPair.first->end())
        slidingFitResultMapPair.first->erase(microFitToDelete);

    const TwoDSlidingFitResultMap::const_iterator macroFitToDelete(slidingFitResultMapPair.second->find(pClusterToRemove));
    if (macroFitToDelete != slidingFitResultMapPair.second->end())
        slidingFitResultMapPair.second->erase(macroFitToDelete);

    const ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pClusterToRemove));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRefinementBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MicroSlidingFitWindow", m_microSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MacroSlidingFitWindow", m_macroSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "StableRegionClusterFraction", m_stableRegionClusterFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MergePointMinCosAngleDeviation", m_mergePointMinCosAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitFractionForHitRemoval", m_minHitFractionForHitRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDistanceFromMainTrack", m_maxDistanceFromMainTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxHitDistanceFromCluster", m_maxHitDistanceFromCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxHitSeparationForConnectedCluster", m_maxHitSeparationForConnectedCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxTrackGaps", m_maxTrackGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LineSegmentLength", m_lineSegmentLength));

    if (m_lineSegmentLength < std::numeric_limits<float>::epsilon())
    {
        std::cout << "TrackInEMShowerAlgorithm: Line segment length must be positive and nonzero" << std::endl;
        throw STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitWidthMode", m_hitWidthMode));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRefinementBaseAlgorithm::SortByDistanceAlongLine::operator()(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs) const
{
    const CartesianVector lhsHitPosition(
        m_hitWidthMode ? LArHitWidthHelper::GetClosestPointToLine2D(m_startPoint, m_lineDirection, pLhs) : pLhs->GetPositionVector());

    const CartesianVector rhsHitPosition(
        m_hitWidthMode ? LArHitWidthHelper::GetClosestPointToLine2D(m_startPoint, m_lineDirection, pRhs) : pRhs->GetPositionVector());

    const float lhsLongitudinal = m_lineDirection.GetDotProduct(lhsHitPosition - m_startPoint);
    const float rhsLongitudinal = m_lineDirection.GetDotProduct(rhsHitPosition - m_startPoint);

    if (std::fabs(lhsLongitudinal - rhsLongitudinal) > std::numeric_limits<float>::epsilon())
    {
        // ATTN: Order from closest to furthest away
        return (lhsLongitudinal < rhsLongitudinal);
    }

    // ATTN: To deal with tiebreaks
    return LArClusterHelper::SortHitsByPulseHeight(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

typedef bool (*SortFunction)(const Cluster *, const Cluster *);

template void TrackRefinementBaseAlgorithm::UpdateContainers<SortFunction>(
    const ClusterList &, const ClusterList &, const SortFunction, ClusterVector &, SlidingFitResultMapPair &) const;
template void TrackRefinementBaseAlgorithm::InitialiseContainers<SortFunction>(
    const ClusterList *, const SortFunction, ClusterVector &, SlidingFitResultMapPair &) const;

} // namespace lar_content
