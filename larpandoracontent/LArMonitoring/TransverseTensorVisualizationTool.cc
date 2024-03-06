/**
 *  @file   larpandoracontent/LArMonitoring/TransverseTensorVisualizationTool.cc
 *
 *  @brief  Implementation of the transverse tensor visualization tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TransverseTensorVisualizationTool.h"

using namespace pandora;

namespace lar_content
{

TransverseTensorVisualizationTool::TransverseTensorVisualizationTool() :
    m_minClusterConnections(1), m_ignoreUnavailableClusters(true), m_showEachIndividualElement(false), m_showContext(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseTensorVisualizationTool::Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterSet usedKeyClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (m_ignoreUnavailableClusters && !pKeyCluster->IsAvailable())
            continue;

        if (usedKeyClusters.count(pKeyCluster))
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, m_ignoreUnavailableClusters, elementList, nU, nV, nW);

        if ((nU < m_minClusterConnections) && (nV < m_minClusterConnections) && (nW < m_minClusterConnections))
            continue;

        if (nU * nV * nW == 0)
            continue;

        int counter(0);
        ClusterList allClusterListU, allClusterListV, allClusterListW;
        std::cout << " Connections: nU " << nU << ", nV " << nV << ", nW " << nW << ", nElements " << elementList.size() << std::endl;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (allClusterListU.end() == std::find(allClusterListU.begin(), allClusterListU.end(), eIter->GetClusterU()))
                allClusterListU.push_back(eIter->GetClusterU());
            if (allClusterListV.end() == std::find(allClusterListV.begin(), allClusterListV.end(), eIter->GetClusterV()))
                allClusterListV.push_back(eIter->GetClusterV());
            if (allClusterListW.end() == std::find(allClusterListW.begin(), allClusterListW.end(), eIter->GetClusterW()))
                allClusterListW.push_back(eIter->GetClusterW());
            usedKeyClusters.insert(eIter->GetClusterU());

            std::cout << " Element " << counter++ << ": MatchedFraction " << eIter->GetOverlapResult().GetMatchedFraction()
                      << ", MatchedSamplingPoints " << eIter->GetOverlapResult().GetNMatchedSamplingPoints() << ", xSpanU "
                      << eIter->GetOverlapResult().GetXOverlap().GetXSpanU() << ", xSpanV " << eIter->GetOverlapResult().GetXOverlap().GetXSpanV()
                      << ", xSpanW " << eIter->GetOverlapResult().GetXOverlap().GetXSpanW() << ", xOverlapSpan "
                      << eIter->GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;

            if (m_showEachIndividualElement)
            {
                const ClusterList clusterListU(1, eIter->GetClusterU()), clusterListV(1, eIter->GetClusterV()),
                    clusterListW(1, eIter->GetClusterW());
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListU, "UCluster", RED));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListV, "VCluster", GREEN));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListW, "WCluster", BLUE));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
        }

        std::cout << " All Connected Clusters " << std::endl;
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterListU, "AllUClusters", RED));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterListV, "AllVClusters", GREEN));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterListW, "AllWClusters", BLUE));

        if (m_showContext)
        {
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_U)), "InputClusterListU", GRAY));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_V)), "InputClusterListV", GRAY));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_W)), "InputClusterListW", GRAY));
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseTensorVisualizationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterConnections", m_minClusterConnections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShowEachIndividualElement", m_showEachIndividualElement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowContext", m_showContext));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
