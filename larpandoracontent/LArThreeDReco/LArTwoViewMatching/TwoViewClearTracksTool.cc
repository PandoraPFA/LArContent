/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewClearTracksTool.cc
 *
 *  @brief  Implementation of the two view clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewClearTracksTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewClearTracksTool::TwoViewClearTracksTool() :
    m_minXOverlapFraction(0.1f),
    m_minLocallyMatchedFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewClearTracksTool::Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
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

void TwoViewClearTracksTool::CreateThreeDParticles(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList,
    bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (MatrixType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        if (iter->GetOverlapResult().GetTwoViewXOverlap().GetXOverlapFraction0() - m_minXOverlapFraction < std::numeric_limits<float>::epsilon())
            continue;
        if (iter->GetOverlapResult().GetTwoViewXOverlap().GetXOverlapFraction1() - m_minXOverlapFraction < std::numeric_limits<float>::epsilon())
            continue;

        if (iter->GetOverlapResult().GetLocallyMatchedFraction() - m_minLocallyMatchedFraction < std::numeric_limits<float>::epsilon())
            continue;

        // TODO Add real logic here

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(iter->GetCluster1());
        protoParticle.m_clusterList.push_back(iter->GetCluster2());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewClearTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
