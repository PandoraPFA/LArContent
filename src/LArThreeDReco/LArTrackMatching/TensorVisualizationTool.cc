/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/TensorVisualizationTool.cc
 * 
 *  @brief  Implementation of the tensor visualization tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/TensorVisualizationTool.h"

using namespace pandora;

namespace lar
{

StatusCode TensorVisualizationTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
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

        int counter(0);
        ClusterList clusterListU, clusterListV, clusterListW;
        std::cout << " Connections: nU " << nU << ", nV " << nV << ", nW " << nW << ", nElements " << elementList.size() << std::endl;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            clusterListU.insert(eIter->GetClusterU());
            clusterListV.insert(eIter->GetClusterV());
            clusterListW.insert(eIter->GetClusterW());
            usedUClusters.insert(eIter->GetClusterU());

            std::cout << " Element " << counter++ << ": MatchedFraction " << eIter->GetOverlapResult().GetMatchedFraction()
                      << ", MatchedSamplingPoints " << eIter->GetOverlapResult().GetNMatchedSamplingPoints()
                      << ", xSpanU " << eIter->GetOverlapResult().GetXOverlap().GetXSpanU()
                      << ", xSpanV " << eIter->GetOverlapResult().GetXOverlap().GetXSpanV()
                      << ", xSpanW " << eIter->GetOverlapResult().GetXOverlap().GetXSpanW()
                      << ", xOverlapSpan " << eIter->GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;
        }

        PANDORA_MONITORING_API(VisualizeClusters(&clusterListU, "UClusters", RED));
        PANDORA_MONITORING_API(VisualizeClusters(&clusterListV, "VClusters", GREEN));
        PANDORA_MONITORING_API(VisualizeClusters(&clusterListW, "WClusters", BLUE));
        PANDORA_MONITORING_API(ViewEvent());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TensorVisualizationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterConnections = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterConnections", m_minClusterConnections));

    m_ignoreUnavailableClusters = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
