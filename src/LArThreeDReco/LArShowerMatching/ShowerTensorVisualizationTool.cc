/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.cc
 * 
 *  @brief  Implementation of the shower tensor visualization tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.h"

using namespace pandora;

namespace lar
{

bool ShowerTensorVisualizationTool::Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ClusterList usedUClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (m_ignoreUnavailableClusters && !iterU->first->IsAvailable())
            continue;

        if (usedUClusters.count(iterU->first))
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, m_ignoreUnavailableClusters, elementList, nU, nV, nW);

        if ((nU < m_minClusterConnections) && (nV < m_minClusterConnections) && (nW < m_minClusterConnections))
            continue;

        if (nU * nV * nW == 0)
            continue;

        int counter(0);
        ClusterList allClusterListU, allClusterListV, allClusterListW;
        std::cout << " Connections: nU " << nU << ", nV " << nV << ", nW " << nW << ", nElements " << elementList.size() << std::endl;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            allClusterListU.insert(eIter->GetClusterU());
            allClusterListV.insert(eIter->GetClusterV());
            allClusterListW.insert(eIter->GetClusterW());
            usedUClusters.insert(eIter->GetClusterU());

            std::cout << " Element " << counter++ << ": MatchedFraction " << eIter->GetOverlapResult().GetMatchedFraction()
                      << ", MatchedSamplingPoints " << eIter->GetOverlapResult().GetNMatchedSamplingPoints()
                      << ", xSpanU " << eIter->GetOverlapResult().GetXOverlap().GetXSpanU()
                      << ", xSpanV " << eIter->GetOverlapResult().GetXOverlap().GetXSpanV()
                      << ", xSpanW " << eIter->GetOverlapResult().GetXOverlap().GetXSpanW()
                      << ", xOverlapSpan " << eIter->GetOverlapResult().GetXOverlap().GetXOverlapSpan()
                      << ", Availability (" << eIter->GetClusterU()->IsAvailable() << eIter->GetClusterV()->IsAvailable() << eIter->GetClusterW()->IsAvailable() << ") "
                      << std::endl;

            if (m_showEachIndividualElement)
            {
                ClusterList clusterListU, clusterListV, clusterListW;
                clusterListU.insert(eIter->GetClusterU());
                clusterListV.insert(eIter->GetClusterV());
                clusterListW.insert(eIter->GetClusterW());

                PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
                PANDORA_MONITORING_API(VisualizeClusters(&clusterListU, "UCluster", RED));
                PANDORA_MONITORING_API(VisualizeClusters(&clusterListV, "VCluster", GREEN));
                PANDORA_MONITORING_API(VisualizeClusters(&clusterListW, "WCluster", BLUE));
                PANDORA_MONITORING_API(ViewEvent());
            }
        }

        std::cout << " All Connected Clusters " << std::endl;
        PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
        PANDORA_MONITORING_API(VisualizeClusters(&allClusterListU, "AllUClusters", RED));
        PANDORA_MONITORING_API(VisualizeClusters(&allClusterListV, "AllVClusters", GREEN));
        PANDORA_MONITORING_API(VisualizeClusters(&allClusterListW, "AllWClusters", BLUE));

        if (m_showContext)
        {
            PANDORA_MONITORING_API(VisualizeClusters(&(pAlgorithm->GetInputClusterListU()), "InputClusterListU", GRAY));
            PANDORA_MONITORING_API(VisualizeClusters(&(pAlgorithm->GetInputClusterListV()), "InputClusterListV", GRAY));
            PANDORA_MONITORING_API(VisualizeClusters(&(pAlgorithm->GetInputClusterListW()), "InputClusterListW", GRAY));
        }

        PANDORA_MONITORING_API(ViewEvent());
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerTensorVisualizationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterConnections = 1;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterConnections", m_minClusterConnections));

    m_ignoreUnavailableClusters = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    m_showEachIndividualElement = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowEachIndividualElement", m_showEachIndividualElement));

    m_showContext = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowContext", m_showContext));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
