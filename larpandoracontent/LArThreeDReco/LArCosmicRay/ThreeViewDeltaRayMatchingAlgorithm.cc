/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeViewDeltaRayMatchingAlgorithm::ThreeViewDeltaRayMatchingAlgorithm() :
    m_minClusterCaloHits(5),
    m_nMaxTensorToolRepeats(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewDeltaRayMatchingAlgorithm::DoesClusterPassTensorThreshold(const Cluster *const pCluster) const
{
    return (pCluster->GetNCaloHits() >= m_minClusterCaloHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    DeltaRayOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult));

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV,
    const Cluster *const pClusterW, DeltaRayOverlapResult &overlapResult) const
{
    float chiSquaredSum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    XOverlap xOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

    StatusCode statusCode(
        this->PerformThreeViewMatching(pClusterU, pClusterV, pClusterW, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xOverlapObject));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    PfoList commonMuonPfoList;
    this->FindCommonMuonParents(pClusterU, pClusterV, pClusterW, commonMuonPfoList);

    if (commonMuonPfoList.empty())
        return STATUS_CODE_NOT_FOUND;

    overlapResult = DeltaRayOverlapResult(nMatchedSamplingPoints, nSamplingPoints, chiSquaredSum, xOverlapObject, commonMuonPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::FindCommonMuonParents(
    const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW, PfoList &commonMuonPfoList) const
{
    ClusterList consideredClustersU, consideredClustersV, consideredClustersW;
    PfoList nearbyMuonPfosU, nearbyMuonPfosV, nearbyMuonPfosW;

    this->GetNearbyMuonPfos(pClusterU, consideredClustersU, nearbyMuonPfosU);

    if (nearbyMuonPfosU.empty())
        return;

    this->GetNearbyMuonPfos(pClusterV, consideredClustersV, nearbyMuonPfosV);

    if (nearbyMuonPfosV.empty())
        return;

    this->GetNearbyMuonPfos(pClusterW, consideredClustersW, nearbyMuonPfosW);

    if (nearbyMuonPfosW.empty())
        return;

    for (const ParticleFlowObject *const pNearbyMuonU : nearbyMuonPfosU)
    {
        for (const ParticleFlowObject *const pNearbyMuonV : nearbyMuonPfosV)
        {
            if (pNearbyMuonV != pNearbyMuonU)
                continue;

            for (const ParticleFlowObject *const pNearbyMuonW : nearbyMuonPfosW)
            {
                if (pNearbyMuonW == pNearbyMuonV)
                    commonMuonPfoList.emplace_back(pNearbyMuonU);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    // Apply tools sequentially restarting if a change is made and ending if the tools finish or the restart limit is reached
    unsigned int repeatCounter(0);

    for (auto toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end();)
    {
        DeltaRayTensorTool *const pTool(*toolIter);
        const bool repeatTools(pTool->Run(this, this->GetMatchingControl().GetOverlapTensor()));

        toolIter = repeatTools ? m_algorithmToolVector.begin() : toolIter + 1;
        repeatCounter = repeatTools ? repeatCounter + 1 : repeatCounter;

        if (repeatCounter > m_nMaxTensorToolRepeats)
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "DeltaRayTools", algorithmToolVector));

    for (auto algorithmTool : algorithmToolVector)
    {
        DeltaRayTensorTool *const pDeltaRayTensorTool(dynamic_cast<DeltaRayTensorTool *>(algorithmTool));

        if (!pDeltaRayTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayTensorTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusterRebuilding", m_reclusteringAlgorithmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
