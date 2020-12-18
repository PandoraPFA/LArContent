/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ClearDeltaRayTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ClearDeltaRayTool.h"

using namespace pandora;

namespace lar_content
{

ClearDeltaRayTool::ClearDeltaRayTool() :
    m_minMatchedFraction(0.9f),
    m_minXOverlapFraction(0.7f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

    bool ClearDeltaRayTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
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

void ClearDeltaRayTool::CreateThreeDParticles(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
    bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        if (iter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        const XOverlap &xOverlap(iter->GetOverlapResult().GetXOverlap());

        if ((xOverlap.GetXSpanU() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanU() < m_minXOverlapFraction))
            continue;

        if ((xOverlap.GetXSpanV() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanV() < m_minXOverlapFraction))
            continue;

        if ((xOverlap.GetXSpanW() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanW() < m_minXOverlapFraction))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(iter->GetClusterU());
        protoParticle.m_clusterList.push_back(iter->GetClusterV());
        protoParticle.m_clusterList.push_back(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreatePfos(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearDeltaRayTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
