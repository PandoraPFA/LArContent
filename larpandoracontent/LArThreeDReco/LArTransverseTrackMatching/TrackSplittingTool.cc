/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.cc
 *
 *  @brief  Implementation of the long tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h"

using namespace pandora;

namespace lar_content
{

TrackSplittingTool::TrackSplittingTool() :
    m_minMatchedFraction(0.75f),
    m_minMatchedSamplingPoints(10),
    m_minXOverlapFraction(0.75f),
    m_minMatchedSamplingPointRatio(2),
    m_maxShortDeltaXFraction(0.2f),
    m_maxAbsoluteShortDeltaX(5.f),
    m_minLongDeltaXFraction(0.2f),
    m_minAbsoluteLongDeltaX(1.f),
    m_minSplitToVertexProjection(1.f),
    m_maxSplitVsFitPositionDistance(1.5f),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackSplittingTool::Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    SplitPositionMap splitPositionMap;
    this->FindTracks(pAlgorithm, overlapTensor, splitPositionMap);

    const bool splitsMade(pAlgorithm->MakeClusterSplits(splitPositionMap));
    return splitsMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackSplittingTool::FindTracks(
    ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType &overlapTensor, SplitPositionMap &splitPositionMap) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectElements(elementList, usedClusters, iteratorList);

        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            if (!this->PassesChecks(pAlgorithm, *(*iIter), true, usedClusters, splitPositionMap) &&
                !this->PassesChecks(pAlgorithm, *(*iIter), false, usedClusters, splitPositionMap))
            {
                continue;
            }

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackSplittingTool::SelectElements(const TensorType::ElementList &elementList, const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
            continue;

        const XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());
        const float longSpan(std::max(xOverlap.GetXSpanU(), std::max(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float shortSpan1(std::min(xOverlap.GetXSpanU(), std::min(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float shortSpan2(
            ((xOverlap.GetXSpanU() > shortSpan1) && (xOverlap.GetXSpanU() < longSpan))
                ? xOverlap.GetXSpanU()
                : ((xOverlap.GetXSpanV() > shortSpan1) && (xOverlap.GetXSpanV() < longSpan)) ? xOverlap.GetXSpanV() : xOverlap.GetXSpanW());

        if ((shortSpan1 < std::numeric_limits<float>::epsilon()) || (longSpan < std::numeric_limits<float>::epsilon()))
            continue;

        if ((shortSpan1 / xOverlap.GetXOverlapSpan()) < m_minXOverlapFraction)
            continue;

        if ((shortSpan1 / shortSpan2) < m_minXOverlapFraction)
            continue;

        if ((shortSpan1 / longSpan) > m_minXOverlapFraction)
            continue;

        iteratorList.push_back(eIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackSplittingTool::PassesChecks(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const bool isMinX, ClusterSet &usedClusters, SplitPositionMap &splitPositionMap) const
{
    try
    {
        const Particle particle(element);

        if (usedClusters.count(particle.m_pLongCluster) || usedClusters.count(particle.m_pCluster1) || usedClusters.count(particle.m_pCluster2))
            return false;

        const float longXSpan(particle.m_longMaxX - particle.m_longMinX);

        if (longXSpan < std::numeric_limits<float>::epsilon())
            return false;

        float splitX(isMinX ? (0.5f * (particle.m_short1MinX + particle.m_short2MinX)) : (0.5f * (particle.m_short1MaxX + particle.m_short2MaxX)));
        const float shortDeltaX(
            isMinX ? std::fabs(particle.m_short1MinX - particle.m_short2MinX) : std::fabs(particle.m_short1MaxX - particle.m_short2MaxX));
        const float longDeltaX(isMinX ? (splitX - particle.m_longMinX) : (particle.m_longMaxX - splitX));

        if (((shortDeltaX / longXSpan) > m_maxShortDeltaXFraction) || (shortDeltaX > m_maxAbsoluteShortDeltaX) ||
            ((longDeltaX / longXSpan) < m_minLongDeltaXFraction) || (longDeltaX < m_minAbsoluteLongDeltaX))
        {
            return false;
        }

        const LArPointingCluster pointingCluster1(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCluster1));
        const LArPointingCluster pointingCluster2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCluster2));
        const LArPointingCluster longPointingCluster(pAlgorithm->GetCachedSlidingFitResult(particle.m_pLongCluster));

        const CartesianVector &minPos1{
            pointingCluster1.GetInnerVertex().GetPosition().GetX() < pointingCluster1.GetOuterVertex().GetPosition().GetX()
                ? pointingCluster1.GetInnerVertex().GetPosition()
                : pointingCluster1.GetOuterVertex().GetPosition()};
        const CartesianVector &minPos2{
            pointingCluster2.GetInnerVertex().GetPosition().GetX() < pointingCluster2.GetOuterVertex().GetPosition().GetX()
                ? pointingCluster2.GetInnerVertex().GetPosition()
                : pointingCluster2.GetOuterVertex().GetPosition()};
        const CartesianVector &maxPos1{
            pointingCluster1.GetInnerVertex().GetPosition().GetX() > pointingCluster1.GetOuterVertex().GetPosition().GetX()
                ? pointingCluster1.GetInnerVertex().GetPosition()
                : pointingCluster1.GetOuterVertex().GetPosition()};
        const CartesianVector &maxPos2{
            pointingCluster2.GetInnerVertex().GetPosition().GetX() > pointingCluster2.GetOuterVertex().GetPosition().GetX()
                ? pointingCluster2.GetInnerVertex().GetPosition()
                : pointingCluster2.GetOuterVertex().GetPosition()};
        CartesianVector position1(isMinX ? minPos1 : maxPos1);
        CartesianVector position2(isMinX ? minPos2 : maxPos2);

        if (!isMinX && (maxPos2.GetX() - maxPos1.GetX() > 1.f))
            this->RefineSplitPosition(pAlgorithm, particle.m_pCluster1, particle.m_pCluster2, position1, position2, splitX);
        else if (!isMinX && (maxPos1.GetX() - maxPos2.GetX() > 1.f))
            this->RefineSplitPosition(pAlgorithm, particle.m_pCluster2, particle.m_pCluster1, position2, position1, splitX);
        else if (isMinX && (minPos2.GetX() - minPos1.GetX() > 1.f))
            this->RefineSplitPosition(pAlgorithm, particle.m_pCluster2, particle.m_pCluster1, position2, position1, splitX);
        else if (isMinX && (minPos1.GetX() - minPos2.GetX() > 1.f))
            this->RefineSplitPosition(pAlgorithm, particle.m_pCluster1, particle.m_pCluster2, position1, position2, splitX);

        CartesianVector splitPosition(0.f, 0.f, 0.f);
        float chiSquared(std::numeric_limits<float>::max());
        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), LArClusterHelper::GetClusterHitType(particle.m_pCluster1),
            LArClusterHelper::GetClusterHitType(particle.m_pCluster2), position1, position2, splitPosition, chiSquared);

        if (m_visualize)
        {
            ClusterList clusterList{particle.m_pLongCluster, particle.m_pCluster1, particle.m_pCluster2};
            PANDORA_MONITORING_API(SetEveDisplayParameters(pAlgorithm->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Clusters", RED));
            PANDORA_MONITORING_API(AddMarkerToVisualization(pAlgorithm->GetPandora(), &position1, "Pos 1", BLACK, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(pAlgorithm->GetPandora(), &position2, "Pos 2", BLACK, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(pAlgorithm->GetPandora(), &splitPosition, "Split Pos", BLACK, 1));
            PANDORA_MONITORING_API(ViewEvent(pAlgorithm->GetPandora()));
        }

        if (!this->CheckSplitPosition(splitPosition, splitX, pAlgorithm->GetCachedSlidingFitResult(particle.m_pLongCluster)))
            return false;

        const CartesianVector splitToInnerVertex(splitPosition - longPointingCluster.GetInnerVertex().GetPosition());
        const CartesianVector outerVertexToSplit(longPointingCluster.GetOuterVertex().GetPosition() - splitPosition);
        const CartesianVector outerToInnerUnitVector(
            (longPointingCluster.GetOuterVertex().GetPosition() - longPointingCluster.GetInnerVertex().GetPosition()).GetUnitVector());

        if ((splitToInnerVertex.GetDotProduct(outerToInnerUnitVector) > m_minSplitToVertexProjection) &&
            (outerVertexToSplit.GetDotProduct(outerToInnerUnitVector) > m_minSplitToVertexProjection))
        {
            splitPositionMap[particle.m_pLongCluster].push_back(splitPosition);
            return true;
        }
    }
    catch (StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackSplittingTool::CheckSplitPosition(const CartesianVector &splitPosition, const float splitX, const TwoDSlidingFitResult &longFitResult) const
{
    CartesianPointVector fitPositionList;
    if (STATUS_CODE_SUCCESS != longFitResult.GetGlobalFitPositionListAtX(splitX, fitPositionList))
        return false;

    for (CartesianPointVector::const_iterator iter = fitPositionList.begin(), iterEnd = fitPositionList.end(); iter != iterEnd; ++iter)
    {
        if ((splitPosition - *iter).GetMagnitude() < m_maxSplitVsFitPositionDistance)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackSplittingTool::RefineSplitPosition(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const Cluster *pFixedCluster,
    const Cluster *pRefineCluster, const CartesianVector &fixedPosition, CartesianVector &refinePosition, float &splitX) const
{
    CaloHitList caloHits;
    pRefineCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
    float minChi2{std::numeric_limits<float>::max()};
    const CaloHit *pBestCaloHit{nullptr};
    for (const CaloHit *pCaloHit : caloHits)
    {
        float chi2{std::numeric_limits<float>::max()};
        CartesianVector dummy(0.f, 0.f, 0.f);
        LArGeometryHelper::MergeTwoPositions3D(pAlgorithm->GetPandora(), LArClusterHelper::GetClusterHitType(pFixedCluster),
            LArClusterHelper::GetClusterHitType(pRefineCluster), fixedPosition, pCaloHit->GetPositionVector(), dummy, chi2);
        if (chi2 < minChi2)
        {
            minChi2 = chi2;
            pBestCaloHit = pCaloHit;
        }
    }
    if (pBestCaloHit)
    {
        const CartesianVector &bestPosition{pBestCaloHit->GetPositionVector()};
        refinePosition.SetValues(bestPosition.GetX(), bestPosition.GetY(), bestPosition.GetZ());
        splitX = bestPosition.GetX();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackSplittingTool::Particle::Particle(const TensorType::Element &element)
{
    const XOverlap &xOverlap(element.GetOverlapResult().GetXOverlap());

    const HitType longHitType =
        ((xOverlap.GetXSpanU() > xOverlap.GetXSpanV()) && (xOverlap.GetXSpanU() > xOverlap.GetXSpanW()))
            ? TPC_VIEW_U
            : ((xOverlap.GetXSpanV() > xOverlap.GetXSpanU()) && (xOverlap.GetXSpanV() > xOverlap.GetXSpanW()))
                  ? TPC_VIEW_V
                  : ((xOverlap.GetXSpanW() > xOverlap.GetXSpanU()) && (xOverlap.GetXSpanW() > xOverlap.GetXSpanV())) ? TPC_VIEW_W : HIT_CUSTOM;

    if (HIT_CUSTOM == longHitType)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_pLongCluster =
        (TPC_VIEW_U == longHitType) ? element.GetClusterU() : (TPC_VIEW_V == longHitType) ? element.GetClusterV() : element.GetClusterW();
    m_pCluster1 = (TPC_VIEW_U == longHitType) ? element.GetClusterV() : element.GetClusterU();
    m_pCluster2 = (TPC_VIEW_W == longHitType) ? element.GetClusterV() : element.GetClusterW();
    m_longMinX = (TPC_VIEW_U == longHitType) ? xOverlap.GetUMinX() : (TPC_VIEW_V == longHitType) ? xOverlap.GetVMinX() : xOverlap.GetWMinX();
    m_longMaxX = (TPC_VIEW_U == longHitType) ? xOverlap.GetUMaxX() : (TPC_VIEW_V == longHitType) ? xOverlap.GetVMaxX() : xOverlap.GetWMaxX();
    m_short1MinX = (TPC_VIEW_U == longHitType) ? xOverlap.GetVMinX() : xOverlap.GetUMinX();
    m_short1MaxX = (TPC_VIEW_U == longHitType) ? xOverlap.GetVMaxX() : xOverlap.GetUMaxX();
    m_short2MinX = (TPC_VIEW_W == longHitType) ? xOverlap.GetVMinX() : xOverlap.GetWMinX();
    m_short2MaxX = (TPC_VIEW_W == longHitType) ? xOverlap.GetVMaxX() : xOverlap.GetWMaxX();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackSplittingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxShortDeltaXFraction", m_maxShortDeltaXFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAbsoluteShortDeltaX", m_maxAbsoluteShortDeltaX));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinLongDeltaXFraction", m_minLongDeltaXFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinAbsoluteLongDeltaX", m_minAbsoluteLongDeltaX));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSplitToVertexProjection", m_minSplitToVertexProjection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxSplitVsFitPositionDistance", m_maxSplitVsFitPositionDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
