/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/MissingTrackTool.cc
 * 
 *  @brief  Implementation of the missing track tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/MissingTrackTool.h"

using namespace pandora;

namespace lar
{

bool MissingTrackTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindMissingTracks(pAlgorithm, overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackTool::FindMissingTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, false, elementList, nU, nV, nW);

        unsigned int counter(0);
        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            unsigned int nAvailable(0);

            if (eIter->GetClusterU()->IsAvailable())
                ++nAvailable;

            if (eIter->GetClusterV()->IsAvailable())
                ++nAvailable;

            if (eIter->GetClusterW()->IsAvailable())
                ++nAvailable;

            if (2 != nAvailable)
                continue;

            // TODO

//std::cout << " Element " << counter++ << ": MatchedFraction " << eIter->GetOverlapResult().GetMatchedFraction()
//<< ", MatchedSamplingPoints " << eIter->GetOverlapResult().GetNMatchedSamplingPoints()
//<< ", xSpanU " << eIter->GetOverlapResult().GetXOverlap().GetXSpanU()
//<< ", xSpanV " << eIter->GetOverlapResult().GetXOverlap().GetXSpanV()
//<< ", xSpanW " << eIter->GetOverlapResult().GetXOverlap().GetXSpanW()
//<< ", xOverlapSpan " << eIter->GetOverlapResult().GetXOverlap().GetXOverlapSpan()
//<< " AU " << eIter->GetClusterU()->IsAvailable()
//<< " AV " << eIter->GetClusterV()->IsAvailable()
//<< " AW " << eIter->GetClusterW()->IsAvailable() << std::endl;
//
//ClusterList clusterListU, clusterListV, clusterListW;
//clusterListU.insert(eIter->GetClusterU());
//clusterListV.insert(eIter->GetClusterV());
//clusterListW.insert(eIter->GetClusterW());
//
//PANDORA_MONITORING_API(SetEveDisplayParameters(false, DETECTOR_VIEW_XZ));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListU, "UCluster", RED));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListV, "VCluster", GREEN));
//PANDORA_MONITORING_API(VisualizeClusters(&clusterListW, "WCluster", BLUE));
//PANDORA_MONITORING_API(VisualizeClusters(&(pAlgorithm->GetInputClusterListU()), "InputClusterListU", GRAY));
//PANDORA_MONITORING_API(VisualizeClusters(&(pAlgorithm->GetInputClusterListV()), "InputClusterListV", GRAY));
//PANDORA_MONITORING_API(VisualizeClusters(&(pAlgorithm->GetInputClusterListW()), "InputClusterListW", GRAY));
//PANDORA_MONITORING_API(ViewEvent());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTrackTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
