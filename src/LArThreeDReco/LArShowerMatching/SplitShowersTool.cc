/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerMatching/SplitShowersTool.cc
 * 
 *  @brief  Implementation of the split showers tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArShowerMatching/SplitShowersTool.h"

using namespace pandora;

namespace lar
{

bool SplitShowersTool::Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ClusterMergeMap clusterMergeMap;
    this->FindSplitShowers(pAlgorithm, overlapTensor, clusterMergeMap);

    return this->ApplyChanges(pAlgorithm, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::FindSplitShowers(ThreeDShowersAlgorithm *pAlgorithm, const TensorType &overlapTensor, ClusterMergeMap &clusterMergeMap) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        if (nU * nV * nW < 2)
            continue;

        std::sort(elementList.begin(), elementList.end(), ThreeDShowersAlgorithm::SortByNMatchedSamplingPoints);

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (!this->PassesElementCuts(eIter, usedClusters))
                continue;

            IteratorList iteratorList;
            this->SelectTensorElements(eIter, elementList, usedClusters, iteratorList);

            if (iteratorList.size() < 2)
                continue;

            this->FindShowerMerges(pAlgorithm, iteratorList, clusterMergeMap);

            for (ClusterMergeMap::const_iterator cIter = clusterMergeMap.begin(), cIterEnd = clusterMergeMap.end(); cIter != cIterEnd; ++cIter)
            {
                usedClusters.insert(cIter->first);
                usedClusters.insert(cIter->second.begin(), cIter->second.end());
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter, const ClusterList &usedClusters) const
{
    if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
        return false;

    if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
        return false;

    if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
    const ClusterList &usedClusters, IteratorList &iteratorList) const
{
    iteratorList.push_back(eIter);

    for (TensorType::ElementList::const_iterator eIter2 = elementList.begin(); eIter2 != elementList.end(); ++eIter2)
    {
        if (eIter == eIter2)
            continue;

        if (!this->PassesElementCuts(eIter2, usedClusters))
            continue;

        for (IteratorList::const_iterator iIter = iteratorList.begin(); iIter != iteratorList.end(); ++iIter)
        {
            if ((*iIter) == eIter2)
                continue;

            unsigned int nMatchedClusters(0);

            if ((*iIter)->GetClusterU() == eIter2->GetClusterU())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterV() == eIter2->GetClusterV())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterW() == eIter2->GetClusterW())
                ++nMatchedClusters;

            if (m_nCommonClusters == nMatchedClusters)
            {
                iteratorList.push_back(eIter2);
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::FindShowerMerges(ThreeDShowersAlgorithm */*pAlgorithm*/, const IteratorList &iteratorList, ClusterMergeMap &/*clusterMergeMap*/) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2)
        {
            if (iIter1 == iIter2)
                continue;

            try
            {
                

//PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
//ClusterList clusterListU1, clusterListV1, clusterListW1;
//clusterListU1.insert((*iIter1)->GetClusterU());
//clusterListV1.insert((*iIter1)->GetClusterV());
//clusterListW1.insert((*iIter1)->GetClusterW());
//std::cout << "SPLITSHOWERSTOOL: CANDIDATE1 MatchedFraction " << (*iIter1)->GetOverlapResult().GetMatchedFraction()
//<< ", MatchedSamplingPoints " << (*iIter1)->GetOverlapResult().GetNMatchedSamplingPoints()
//<< ", xSpanU " << (*iIter1)->GetOverlapResult().GetXOverlap().GetXSpanU()
//<< ", xSpanV " << (*iIter1)->GetOverlapResult().GetXOverlap().GetXSpanV()
//<< ", xSpanW " << (*iIter1)->GetOverlapResult().GetXOverlap().GetXSpanW()
//<< ", xOverlapSpan " << (*iIter1)->GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListU1, "clusterListU1", RED));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListV1, "clusterListV1", GREEN));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListW1, "clusterListW1", BLUE));
//ClusterList clusterListU2, clusterListV2, clusterListW2;
//clusterListU2.insert((*iIter2)->GetClusterU());
//clusterListV2.insert((*iIter2)->GetClusterV());
//clusterListW2.insert((*iIter2)->GetClusterW());
//std::cout << "SPLITSHOWERSTOOL: CANDIDATE1 MatchedFraction " << (*iIter2)->GetOverlapResult().GetMatchedFraction()
//<< ", MatchedSamplingPoints " << (*iIter2)->GetOverlapResult().GetNMatchedSamplingPoints()
//<< ", xSpanU " << (*iIter2)->GetOverlapResult().GetXOverlap().GetXSpanU()
//<< ", xSpanV " << (*iIter2)->GetOverlapResult().GetXOverlap().GetXSpanV()
//<< ", xSpanW " << (*iIter2)->GetOverlapResult().GetXOverlap().GetXSpanW()
//<< ", xOverlapSpan " << (*iIter2)->GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListU2, "clusterListU2", RED));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListV2, "clusterListV2", GREEN));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListW2, "clusterListW2", BLUE));
//PANDORA_MONITORING_API(ViewEvent());
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::ApplyChanges(ThreeDShowersAlgorithm *pAlgorithm, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterMergeMap consolidatedMergeMap;

    for (ClusterMergeMap::const_iterator cIter = clusterMergeMap.begin(), cIterEnd = clusterMergeMap.end(); cIter != cIterEnd; ++cIter)
    {
        const ClusterList &daughterClusters(cIter->second);

        for (ClusterList::const_iterator dIter = daughterClusters.begin(), dIterEnd = daughterClusters.end(); dIter != dIterEnd; ++dIter)
        {
            if (consolidatedMergeMap.count(*dIter))
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        ClusterList &targetClusterList(consolidatedMergeMap[cIter->first]);
        targetClusterList.insert(daughterClusters.begin(), daughterClusters.end());
    }

    return pAlgorithm->MakeClusterMerges(consolidatedMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SplitShowersTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_nCommonClusters = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NCommonClusters", m_nCommonClusters));

    m_minMatchedFraction = 0.2f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 40;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
