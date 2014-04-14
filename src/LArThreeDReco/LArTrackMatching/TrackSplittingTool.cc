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

    //this->PerformClusterSplits(pAlgorithm, splitPositionMap); TODO
    const bool splitsMade(!splitPositionMap.empty());

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

        if ((shortSpan2 / shortSpan1) < m_minXOverlapFraction)
            continue;

        if ((shortSpan1 / longSpan) > 0.9f) // TODO
            continue;

        iteratorList.push_back(eIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackSplittingTool::PassesChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element,
    ClusterList &usedClusters, SplitPositionMap &splitPositionMap) const
{
    // TODO Get split positions by looking for close end positions of short clusters and deciding they are suitably far away from ends of long cluster

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackSplittingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minXOverlapFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_minMatchedSamplingPointRatio = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
