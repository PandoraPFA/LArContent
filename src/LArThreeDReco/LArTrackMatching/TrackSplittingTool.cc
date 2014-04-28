/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/TrackSplittingTool.cc
 * 
 *  @brief  Implementation of the long tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTrackMatching/LongTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/TrackSplittingTool.h"

using namespace pandora;

namespace lar
{

bool TrackSplittingTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    SplitPositionMap splitPositionMap;
    this->FindTracks(pAlgorithm, overlapTensor, splitPositionMap);

    const bool splitsMade(pAlgorithm->MakeClusterSplits(splitPositionMap));
    return splitsMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackSplittingTool::FindTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, SplitPositionMap &splitPositionMap) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectElements(elementList, usedClusters, iteratorList);

        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            if (!this->PassesChecks(pAlgorithm, *(*iIter), usedClusters, splitPositionMap))
                continue;

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackSplittingTool::SelectElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
            continue;

        const TransverseOverlapResult::XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());
        const float longSpan(std::max(xOverlap.GetXSpanU(), std::max(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float shortSpan1(std::min(xOverlap.GetXSpanU(), std::min(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float shortSpan2(((xOverlap.GetXSpanU() > shortSpan1) && (xOverlap.GetXSpanU() < longSpan)) ? xOverlap.GetXSpanU() :
            ((xOverlap.GetXSpanV() > shortSpan1) && (xOverlap.GetXSpanV() < longSpan)) ? xOverlap.GetXSpanV() : xOverlap.GetXSpanW());

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

bool TrackSplittingTool::PassesChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element,
    ClusterList &usedClusters, SplitPositionMap &splitPositionMap) const
{
    const Particle particle(element);
    const TwoDSlidingFitResult &longFitResult(pAlgorithm->GetCachedSlidingFitResult(particle.m_pLongCluster));
    const float longXSpan(particle.m_longMaxX - particle.m_longMinX);

    if (longXSpan < std::numeric_limits<float>::epsilon())
        return false;

    bool passesChecks(false);

    const float splitMinX(0.5f * (particle.m_short1MinX + particle.m_short2MinX));
    const float shortDeltaMinX(std::fabs(particle.m_short1MinX - particle.m_short2MinX));
    const float longDeltaMinX(splitMinX - particle.m_longMinX);

    if (((shortDeltaMinX / longXSpan) < m_maxShortDeltaXFraction) && (shortDeltaMinX < m_maxAbsoluteShortDeltaX) &&
        ((longDeltaMinX / longXSpan) > m_minLongDeltaXFraction) && (longDeltaMinX > m_minAbsoluteLongDeltaX))
    {
        CartesianVector splitPosition(0.f, 0.f, 0.f);
        longFitResult.GetGlobalFitPosition(splitMinX, true, splitPosition);
        splitPositionMap[particle.m_pLongCluster].push_back(splitPosition);
        passesChecks = true;
    }

    const float splitMaxX(0.5f * (particle.m_short1MaxX + particle.m_short2MaxX));
    const float shortDeltaMaxX(std::fabs(particle.m_short1MaxX - particle.m_short2MaxX));
    const float longDeltaMaxX(particle.m_longMaxX - splitMaxX);

    if (((shortDeltaMaxX / longXSpan) < m_maxShortDeltaXFraction) && (shortDeltaMaxX < m_maxAbsoluteShortDeltaX) &&
        ((longDeltaMaxX / longXSpan) > m_minLongDeltaXFraction) && (longDeltaMaxX > m_minAbsoluteLongDeltaX))
    {
        CartesianVector splitPosition(0.f, 0.f, 0.f);
        longFitResult.GetGlobalFitPosition(splitMaxX, true, splitPosition);
        splitPositionMap[particle.m_pLongCluster].push_back(splitPosition);
        passesChecks = true;
    }

    return passesChecks;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackSplittingTool::Particle::Particle(const TensorType::Element &element)
{
    const TransverseOverlapResult::XOverlap &xOverlap(element.GetOverlapResult().GetXOverlap());

    const HitType longHitType = ((xOverlap.GetXSpanU() > xOverlap.GetXSpanV()) && (xOverlap.GetXSpanU() > xOverlap.GetXSpanW())) ? TPC_VIEW_U :
        ((xOverlap.GetXSpanV() > xOverlap.GetXSpanU()) && (xOverlap.GetXSpanV() > xOverlap.GetXSpanW())) ? TPC_VIEW_V :
        ((xOverlap.GetXSpanW() > xOverlap.GetXSpanU()) && (xOverlap.GetXSpanW() > xOverlap.GetXSpanV())) ? TPC_VIEW_W : CUSTOM;

    if (CUSTOM == longHitType)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_pLongCluster = (TPC_VIEW_U == longHitType) ? element.GetClusterU() : (TPC_VIEW_V == longHitType) ? element.GetClusterV() : element.GetClusterW();
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
    m_minMatchedFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minXOverlapFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_minMatchedSamplingPointRatio = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    m_maxShortDeltaXFraction = 0.1f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShortDeltaXFraction", m_maxShortDeltaXFraction));

    m_maxAbsoluteShortDeltaX = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAbsoluteShortDeltaX", m_maxAbsoluteShortDeltaX));

    m_minLongDeltaXFraction = 0.2f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongDeltaXFraction", m_minLongDeltaXFraction));

    m_minAbsoluteLongDeltaX = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAbsoluteLongDeltaX", m_minAbsoluteLongDeltaX));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
