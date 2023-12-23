/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.h"

using namespace pandora;

namespace lar_content
{

ClearLongitudinalTracksTool::ClearLongitudinalTracksTool() : m_minMatchedFraction(0.8f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearLongitudinalTracksTool::Run(ThreeViewLongitudinalTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool particlesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearLongitudinalTracksTool::CreateThreeDParticles(
    ThreeViewLongitudinalTracksAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        std::cout << " FRACTION: " << iter->GetOverlapResult().GetMatchedFraction() << std::endl;
        if (iter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(iter->GetClusterU());
        protoParticle.m_clusterList.push_back(iter->GetClusterV());
        protoParticle.m_clusterList.push_back(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearLongitudinalTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
