/**
 *  @file   larpandoracontent/LArMonitoring/TransverseMatrixVisualizationTool.cc
 *
 *  @brief  Implementation of the transverse matrix visualization tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TransverseMatrixVisualizationTool.h"

using namespace pandora;

namespace lar_content
{

TransverseMatrixVisualizationTool::TransverseMatrixVisualizationTool() :
    m_minClusterConnections(1),
    m_ignoreUnavailableClusters(true),
    m_showEachIndividualElement(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseMatrixVisualizationTool::Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterSet usedKeyClusters;
    ClusterVector sortedKeyClusters;
    overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (m_ignoreUnavailableClusters && !pKeyCluster->IsAvailable())
            continue;

        if (usedKeyClusters.count(pKeyCluster))
            continue;

        unsigned int n1(0), n2(0);
        MatrixType::ElementList elementList;
        overlapMatrix.GetConnectedElements(pKeyCluster, m_ignoreUnavailableClusters, elementList, n1, n2);

        if ((n1 < m_minClusterConnections) && (n2 < m_minClusterConnections))
            continue;

        if (n1 * n2 == 0)
            continue;

        int counter(0);
        ClusterList allClusterList1, allClusterList2;
        std::cout << " Connections: n1 " << n1 << ", n2 " << n2 << ", nElements " << elementList.size() << std::endl;

        for (MatrixType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (allClusterList1.end() == std::find(allClusterList1.begin(), allClusterList1.end(), eIter->GetCluster1())) allClusterList1.push_back(eIter->GetCluster1());
            if (allClusterList2.end() == std::find(allClusterList2.begin(), allClusterList2.end(), eIter->GetCluster2())) allClusterList2.push_back(eIter->GetCluster2());
            usedKeyClusters.insert(eIter->GetCluster1());

            std::cout << " Element " << counter++ << ": XOverlap " << eIter->GetOverlapResult() << std::endl;

            if (m_showEachIndividualElement)
            {
                const ClusterList clusterList1(1, eIter->GetCluster1()), clusterList2(1, eIter->GetCluster2());
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList1, "Cluster1", RED));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList2, "Cluster2", GREEN));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
        }

        std::cout << " All Connected Clusters " << std::endl;
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterList1, "AllClusters1", RED));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterList2, "AllClusters2", GREEN));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseMatrixVisualizationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterConnections", m_minClusterConnections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowEachIndividualElement", m_showEachIndividualElement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
