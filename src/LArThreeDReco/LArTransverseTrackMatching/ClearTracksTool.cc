/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.cc
 * 
 *  @brief  Implementation of the clear tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h"

using namespace pandora;

namespace lar_content
{

bool ClearTracksTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
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

void ClearTracksTool::CreateThreeDParticles(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList,
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
        protoParticle.m_clusterListU.insert(iter->GetClusterU());
        protoParticle.m_clusterListV.insert(iter->GetClusterV());
        protoParticle.m_clusterListW.insert(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minXOverlapFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
