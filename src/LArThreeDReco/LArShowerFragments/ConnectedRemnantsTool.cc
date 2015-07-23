/**
 *  @file   LArContent/src/LArThreeDReco/ShowerFragments/ConnectedRemnantsTool.cc
 *
 *  @brief  Implementation of the clear remnants tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h"

using namespace pandora;

namespace lar_content
{

ConnectedRemnantsTool::ConnectedRemnantsTool() :
    m_maxClusterSeparation(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConnectedRemnantsTool::Run(ThreeDRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector; ClusterMergeMap clusterMergeMap;
    this->FindConnectedShowers(overlapTensor, protoParticleVector, clusterMergeMap);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    const bool mergesMade(pAlgorithm->MakeClusterMerges(clusterMergeMap));

    return (particlesMade||mergesMade);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectedRemnantsTool::FindConnectedShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector,
    ClusterMergeMap &clusterMergeMap) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        TensorType::ElementList connectedElements;
        overlapTensor.GetConnectedElements(iterU->first, true, connectedElements);

        ClusterList clusterListU, clusterListV, clusterListW;
        this->GetClusters(connectedElements, usedClusters, clusterListU, clusterListV, clusterListW);

        if (clusterListU.empty() || clusterListV.empty() || clusterListW.empty())
            continue;

        usedClusters.insert(clusterListU.begin(), clusterListU.end());
        usedClusters.insert(clusterListV.begin(), clusterListV.end());
        usedClusters.insert(clusterListW.begin(), clusterListW.end());

        if (!(this->IsConnected(clusterListU) && this->IsConnected(clusterListV) && this->IsConnected(clusterListW)))
            continue;

        const Cluster *const pClusterU = *(clusterListU.begin());
        const Cluster *const pClusterV = *(clusterListV.begin());
        const Cluster *const pClusterW = *(clusterListW.begin());

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(pClusterU);
        protoParticle.m_clusterListV.insert(pClusterV);
        protoParticle.m_clusterListW.insert(pClusterW);
        protoParticleVector.push_back(protoParticle);

        this->FillMergeMap(pClusterU, clusterListU, clusterMergeMap);
        this->FillMergeMap(pClusterV, clusterListV, clusterMergeMap);
        this->FillMergeMap(pClusterW, clusterListW, clusterMergeMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectedRemnantsTool::GetClusters(const TensorType::ElementList &elementList, const ClusterList &usedClusters,
    ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        clusterListU.insert(eIter->GetClusterU());
        clusterListV.insert(eIter->GetClusterV());
        clusterListW.insert(eIter->GetClusterW());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectedRemnantsTool::FillMergeMap(const Cluster *const pFirstCluster, const ClusterList &clusterList, ClusterMergeMap &clusterMergeMap) const
{
    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (!clusterList.count(pFirstCluster))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pSecondCluster = *iter;

        if (pFirstCluster == pSecondCluster)
            continue;

        clusterMergeMap[pFirstCluster].insert(pSecondCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConnectedRemnantsTool::IsConnected(const ClusterList &clusterList) const
{
    for (ClusterList::const_iterator iter1 = clusterList.begin(), iterEnd1 = clusterList.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;

        for (ClusterList::const_iterator iter2 = iter1, iterEnd2 = clusterList.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster *const pCluster2 = *iter2;

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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterSeparation", m_maxClusterSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
