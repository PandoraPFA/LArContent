/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerMatching/SimpleShowersTool.cc
 * 
 *  @brief  Implementation of the clear showers tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArShowerMatching/SimpleShowersTool.h"

using namespace pandora;

namespace lar
{

bool SimpleShowersTool::Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindBestShower(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleShowersTool::FindBestShower(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        TensorType::Element bestElement(NULL, NULL, NULL, ShowerOverlapResult());

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (eIter->GetOverlapResult() > bestElement.GetOverlapResult())
                bestElement = *eIter;
        }

        if (!bestElement.GetOverlapResult().IsInitialized())
            continue;

        if ((NULL == bestElement.GetClusterU()) || (NULL == bestElement.GetClusterV()) || (NULL == bestElement.GetClusterW()))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(bestElement.GetClusterU());
        protoParticle.m_clusterListV.insert(bestElement.GetClusterV());
        protoParticle.m_clusterListW.insert(bestElement.GetClusterW());
        protoParticleVector.push_back(protoParticle);

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
    m_minMatchedFraction = 0.2f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 40;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minXOverlapFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
