/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TwoViewClearDeltaRayTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewClearDeltaRayTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewClearDeltaRayTool::TwoViewClearDeltaRayTool() :
    m_minMatchedFraction(0.9f),
    m_minXOverlapFraction(0.6f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewClearDeltaRayTool::Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool particlesMade(false);

    MatrixType::ElementList elementList;
    overlapMatrix.GetUnambiguousElements(true, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewClearDeltaRayTool::CreateThreeDParticles(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList,
    bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (MatrixType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        if (iter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;
        
        const TwoViewXOverlap &xOverlap(iter->GetOverlapResult().GetXOverlap());

        if ((xOverlap.GetXSpan0() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetTwoViewXOverlapSpan() / xOverlap.GetXSpan0() < m_minXOverlapFraction))
            continue;

        if ((xOverlap.GetXSpan1() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetTwoViewXOverlapSpan() / xOverlap.GetXSpan1() < m_minXOverlapFraction))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(iter->GetCluster1());
        protoParticle.m_clusterList.push_back(iter->GetCluster2());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewClearDeltaRayTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
