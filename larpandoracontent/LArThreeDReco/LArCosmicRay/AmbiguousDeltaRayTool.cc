/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/AmbiguousDeltaRayTool.cc
 *
 *  @brief  Implementation of the ambiguous delta ray tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/AmbiguousDeltaRayTool.h"

using namespace pandora;

namespace lar_content
{

AmbiguousDeltaRayTool::AmbiguousDeltaRayTool() :
    m_maxGoodMatchReducedChiSquared(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AmbiguousDeltaRayTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    m_pParentAlgorithm = pAlgorithm;

    if (PandoraContentApi::GetSettings(*m_pParentAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    this->ExamineConnectedElements(overlapTensor);

    // ATTN: Prevent tensor tool loop running again
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousDeltaRayTool::ExamineConnectedElements(TensorType &overlapTensor) const
{
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    ClusterSet usedClusters;
    ClusterSet usedKeyClusters;
    ProtoParticleVector protoParticleVector;

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (usedKeyClusters.count(pKeyCluster))
            continue;

        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList);

        for (const TensorType::Element &element : elementList)
            usedKeyClusters.insert(element.GetClusterU());

        if (elementList.size() < 2)
            continue;

        this->PickOutGoodMatches(elementList, usedClusters, protoParticleVector);
    }

    if (!protoParticleVector.empty())
        m_pParentAlgorithm->CreatePfos(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousDeltaRayTool::PickOutGoodMatches(
    const TensorType::ElementList &elementList, ClusterSet &usedClusters, ProtoParticleVector &protoParticleVector) const
{
    bool found(false);

    do
    {
        found = false;

        unsigned int highestHitCount(0);
        float bestChiSquared(std::numeric_limits<float>::max());
        const Cluster *pBestClusterU(nullptr), *pBestClusterV(nullptr), *pBestClusterW(nullptr);

        for (const TensorType::Element &element : elementList)
        {
            const Cluster *const pClusterU(element.GetClusterU()), *const pClusterV(element.GetClusterV()), *const pClusterW(element.GetClusterW());

            if (usedClusters.count(pClusterU) || usedClusters.count(pClusterV) || usedClusters.count(pClusterW))
                continue;

            const float chiSquared(element.GetOverlapResult().GetReducedChi2());

            if (chiSquared > m_maxGoodMatchReducedChiSquared)
                continue;

            const unsigned int hitSum(pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits());

            if ((hitSum > highestHitCount) || ((hitSum == highestHitCount) && (chiSquared < bestChiSquared)))
            {
                bestChiSquared = chiSquared;
                highestHitCount = hitSum;
                pBestClusterU = pClusterU;
                pBestClusterV = pClusterV;
                pBestClusterW = pClusterW;
            }
        }

        if (pBestClusterU && pBestClusterV && pBestClusterW)
        {
            found = true;
            usedClusters.insert(pBestClusterU);
            usedClusters.insert(pBestClusterV);
            usedClusters.insert(pBestClusterW);

            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back(pBestClusterU);
            protoParticle.m_clusterList.push_back(pBestClusterV);
            protoParticle.m_clusterList.push_back(pBestClusterW);
            protoParticleVector.push_back(protoParticle);
        }
    } while (found);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AmbiguousDeltaRayTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxGoodMatchReducedChiSquared", m_maxGoodMatchReducedChiSquared));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
