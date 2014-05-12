/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/MatchedEndPointsTool.cc
 *
 *  @brief  Implementation of the matched end point tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/MatchedEndPointsTool.h"

using namespace pandora;

namespace lar
{

bool MatchedEndPointsTool::Run(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindMatchedTracks(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchedEndPointsTool::FindMatchedTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        if (nU * nV * nW == 0)
            continue;


        // --- EVENT DISPLAY [BEGIN] ---
        ClusterList clusterListU, clusterListV, clusterListW;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            clusterListU.insert(eIter->GetClusterU());
            clusterListV.insert(eIter->GetClusterV());
            clusterListW.insert(eIter->GetClusterW());
        }


        PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
        PANDORA_MONITORING_API(VisualizeClusters(&clusterListU, "UCluster", RED));
        PANDORA_MONITORING_API(VisualizeClusters(&clusterListV, "VCluster", GREEN));
        PANDORA_MONITORING_API(VisualizeClusters(&clusterListW, "WCluster", BLUE));
        PANDORA_MONITORING_API(ViewEvent());

        // --- EVENT DISPLAY [END] ---
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchedEndPointsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
