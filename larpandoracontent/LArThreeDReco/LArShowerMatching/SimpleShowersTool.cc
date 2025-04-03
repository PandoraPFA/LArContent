/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/SimpleShowersTool.cc
 *
 *  @brief  Implementation of the clear showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArShowerMatching/SimpleShowersTool.h"

using namespace pandora;

namespace lar_content
{

SimpleShowersTool::SimpleShowersTool() :
    m_minMatchedFraction(0.2f),
    m_minMatchedSamplingPoints(40),
    m_minXOverlapFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SimpleShowersTool::Run(ThreeViewShowersAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindBestShower(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleShowersTool::FindBestShower(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        if (elementList.empty())
            continue;

        TensorType::Element bestElement(elementList.back());

        if (!bestElement.GetOverlapResult().IsInitialized())
            continue;

        if ((NULL == bestElement.GetClusterU()) || (NULL == bestElement.GetClusterV()) || (NULL == bestElement.GetClusterW()))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.emplace_back(bestElement.GetClusterU());
        protoParticle.m_clusterList.emplace_back(bestElement.GetClusterV());
        protoParticle.m_clusterList.emplace_back(bestElement.GetClusterW());
        protoParticleVector.emplace_back(protoParticle);

        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SimpleShowersTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter) const
{
    if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
        return false;

    if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    const XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());

    if ((xOverlap.GetXSpanU() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanU() > m_minXOverlapFraction) &&
        (xOverlap.GetXSpanV() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanV() > m_minXOverlapFraction) &&
        (xOverlap.GetXSpanW() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanW() > m_minXOverlapFraction))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleShowersTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
