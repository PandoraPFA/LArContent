/**
 *  @file   larpandoracontent/LArThreeDReco/ShowerFragments/ConnectedRemnantsTool.cc
 *
 *  @brief  Implementation of the clear remnants tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h"

using namespace pandora;

namespace lar_content
{

ConnectedRemnantsTool::ConnectedRemnantsTool() :
    m_maxClusterSeparation(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConnectedRemnantsTool::Run(ThreeViewRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    ClusterMergeMap clusterMergeMap;
    this->FindConnectedShowers(overlapTensor, protoParticleVector, clusterMergeMap);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    const bool mergesMade(pAlgorithm->MakeClusterMerges(clusterMergeMap));

    return (particlesMade || mergesMade);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectedRemnantsTool::FindConnectedShowers(
    const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        TensorType::ElementList connectedElements;
        overlapTensor.GetConnectedElements(pKeyCluster, true, connectedElements);

        ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;
        this->GetClusters(connectedElements, usedClusters, clusterVectorU, clusterVectorV, clusterVectorW);

        if (clusterVectorU.empty() || clusterVectorV.empty() || clusterVectorW.empty())
            continue;

        usedClusters.insert(clusterVectorU.begin(), clusterVectorU.end());
        usedClusters.insert(clusterVectorV.begin(), clusterVectorV.end());
        usedClusters.insert(clusterVectorW.begin(), clusterVectorW.end());

        if (!(this->IsConnected(clusterVectorU) && this->IsConnected(clusterVectorV) && this->IsConnected(clusterVectorW)))
            continue;

        const Cluster *const pClusterU = clusterVectorU.front();
        const Cluster *const pClusterV = clusterVectorV.front();
        const Cluster *const pClusterW = clusterVectorW.front();

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(pClusterU);
        protoParticle.m_clusterList.push_back(pClusterV);
        protoParticle.m_clusterList.push_back(pClusterW);
        protoParticleVector.push_back(protoParticle);

        this->FillMergeMap(pClusterU, clusterVectorU, clusterMergeMap);
        this->FillMergeMap(pClusterV, clusterVectorV, clusterMergeMap);
        this->FillMergeMap(pClusterW, clusterVectorW, clusterMergeMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectedRemnantsTool::GetClusters(const TensorType::ElementList &elementList, const ClusterSet &usedClusters,
    ClusterVector &clusterVectorU, ClusterVector &clusterVectorV, ClusterVector &clusterVectorW) const
{
    for (const TensorType::Element &element : elementList)
    {
        if (usedClusters.count(element.GetClusterU()) || usedClusters.count(element.GetClusterV()) || usedClusters.count(element.GetClusterW()))
            continue;

        clusterVectorU.push_back(element.GetClusterU());
        clusterVectorV.push_back(element.GetClusterV());
        clusterVectorW.push_back(element.GetClusterW());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectedRemnantsTool::FillMergeMap(const Cluster *const pFirstCluster, const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    if (clusterVector.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (const Cluster *const pSecondCluster : clusterVector)
    {
        if (pFirstCluster == pSecondCluster)
            continue;

        ClusterList &clusterList(clusterMergeMap[pFirstCluster]);

        if (clusterList.end() == std::find(clusterList.begin(), clusterList.end(), pSecondCluster))
            clusterList.push_back(pSecondCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConnectedRemnantsTool::IsConnected(const ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster1 : clusterVector)
    {
        for (const Cluster *const pCluster2 : clusterVector)
        {
            if (pCluster1 == pCluster2)
                continue;

            if (LArClusterHelper::GetClosestDistance(pCluster1, pCluster2) > m_maxClusterSeparation)
                return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConnectedRemnantsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
