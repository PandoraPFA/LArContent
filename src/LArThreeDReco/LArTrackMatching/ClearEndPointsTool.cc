/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ClearEndPointsTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/ClearEndPointsTool.h"

using namespace pandora;

namespace lar
{

bool ClearEndPointsTool::Run(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    bool particlesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearEndPointsTool::CreateThreeDParticles(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList,
    bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
      const float m_reducedChi2Cut = 5.f;  // TODO: Move this to settings

        if (iter->GetOverlapResult().GetReducedChi2() > m_reducedChi2Cut)
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(iter->GetClusterU());
        protoParticle.m_clusterListV.insert(iter->GetClusterV());
        protoParticle.m_clusterListW.insert(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearEndPointsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    // TODO: Fill in settings

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
