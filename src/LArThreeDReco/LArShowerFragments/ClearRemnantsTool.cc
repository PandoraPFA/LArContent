/**
 *  @file   LArContent/src/LArThreeDReco/ShowerFragments/ClearRemnantsTool.cc
 *
 *  @brief  Implementation of the clear remnants tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"

using namespace pandora;

namespace lar_content
{

bool ClearRemnantsTool::Run(ThreeDRemnantsAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    bool particlesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearRemnantsTool::CreateThreeDParticles(ThreeDRemnantsAlgorithm *pAlgorithm, const TensorType::ElementList &elementList,
    bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(iter->GetClusterU());
        protoParticle.m_clusterListV.insert(iter->GetClusterV());
        protoParticle.m_clusterListW.insert(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearRemnantsTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
