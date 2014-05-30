/**
 *  @file   LArContent/src/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.cc
 *
 *  @brief  Implementation of the matched end point tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h"

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

        std::sort(elementList.begin(), elementList.end(), ThreeDLongitudinalTracksAlgorithm::SortByChiSquared);

        for (TensorType::ElementList::const_iterator iter = elementList.begin(); iter != elementList.end(); ++iter)
        {
            if (usedClusters.count(iter->GetClusterU()) || usedClusters.count(iter->GetClusterV()) || usedClusters.count(iter->GetClusterW()))
                continue;

            if (iter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
                continue;

            if (std::max(iter->GetOverlapResult().GetInnerChi2(),iter->GetOverlapResult().GetOuterChi2()) > m_maxEndPointChi2)
                continue;

// --- EVENT DISPLAY [BEGIN] ---
// ClusterList clusterListU, clusterListV, clusterListW;
// clusterListU.insert(iter->GetClusterU());
// clusterListV.insert(iter->GetClusterV());
// clusterListW.insert(iter->GetClusterW());
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
// PANDORA_MONITORING_API(VisualizeClusters(&clusterListU, "ClusterU", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&clusterListV, "ClusterV", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&clusterListW, "ClusterW", BLUE));
// PANDORA_MONITORING_API(ViewEvent());
// --- EVENT DISPLAY [END] ---

            ProtoParticle protoParticle;
            protoParticle.m_clusterListU.insert(iter->GetClusterU());
            protoParticle.m_clusterListV.insert(iter->GetClusterV());
            protoParticle.m_clusterListW.insert(iter->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert(iter->GetClusterU());
            usedClusters.insert(iter->GetClusterV());
            usedClusters.insert(iter->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchedEndPointsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_maxEndPointChi2 = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxEndPointChi2", m_maxEndPointChi2));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
