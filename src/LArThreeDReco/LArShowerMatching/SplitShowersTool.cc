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

    ProtoParticleVector protoParticleVector;
    this->FindSplitShowers(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::FindSplitShowers(const TensorType &overlapTensor, ProtoParticleVector &/*protoParticleVector*/) const
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
ClusterList allClusterListU, allClusterListV, allClusterListW;
        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (!this->PassesElementCuts(eIter, usedClusters))
                continue;

            IteratorList iteratorList;
            this->SelectTensorElements(eIter, elementList, usedClusters, iteratorList);

            if (iteratorList.size() < 2)
                continue;

            for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
            {
allClusterListU.insert((*iIter)->GetClusterU());
allClusterListV.insert((*iIter)->GetClusterV());
allClusterListW.insert((*iIter)->GetClusterW());
std::cout << " Selected Element: MatchedFraction " << (*iIter)->GetOverlapResult().GetMatchedFraction()
<< ", MatchedSamplingPoints " << (*iIter)->GetOverlapResult().GetNMatchedSamplingPoints()
<< ", xSpanU " << (*iIter)->GetOverlapResult().GetXOverlap().GetXSpanU()
<< ", xSpanV " << (*iIter)->GetOverlapResult().GetXOverlap().GetXSpanV()
<< ", xSpanW " << (*iIter)->GetOverlapResult().GetXOverlap().GetXSpanW()
<< ", xOverlapSpan " << (*iIter)->GetOverlapResult().GetXOverlap().GetXOverlapSpan()
<< ", Availability (" << (*iIter)->GetClusterU()->IsAvailable() << (*iIter)->GetClusterV()->IsAvailable() << (*iIter)->GetClusterW()->IsAvailable() << ") "
<< std::endl;
                usedClusters.insert((*iIter)->GetClusterU());
                usedClusters.insert((*iIter)->GetClusterV());
                usedClusters.insert((*iIter)->GetClusterW());
            }
        }
if (!allClusterListU.empty() || !allClusterListV.empty() || !allClusterListW.empty())
{
std::cout << " All Connected Clusters " << std::endl;
PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
PANDORA_MONITORING_API(VisualizeClusters(&allClusterListU, "AllUClusters", RED));
PANDORA_MONITORING_API(VisualizeClusters(&allClusterListV, "AllVClusters", GREEN));
PANDORA_MONITORING_API(VisualizeClusters(&allClusterListW, "AllWClusters", BLUE));
PANDORA_MONITORING_API(ViewEvent());
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
