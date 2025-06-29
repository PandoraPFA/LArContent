/**
 *  @file   larpandoracontent/LArThreeDReco/ShowerFragments/ClearRemnantsTool.cc
 *
 *  @brief  Implementation of the clear remnants tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"

using namespace pandora;

namespace lar_content
{

bool ClearRemnantsTool::Run(ThreeViewRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
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

void ClearRemnantsTool::CreateThreeDParticles(ThreeViewRemnantsAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        ProtoParticle protoParticle;
        protoParticle.m_clusterList.emplace_back(iter->GetClusterU());
        protoParticle.m_clusterList.emplace_back(iter->GetClusterV());
        protoParticle.m_clusterList.emplace_back(iter->GetClusterW());
        protoParticleVector.emplace_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearRemnantsTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
