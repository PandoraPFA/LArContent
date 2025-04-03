/**
 *  @file   larpandoracontent/LArThreeDReco/ShowerFragments/MopUpRemnantsTool.cc
 *
 *  @brief  Implementation of the mop-up remnants tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/MopUpRemnantsTool.h"

using namespace pandora;

namespace lar_content
{

MopUpRemnantsTool::MopUpRemnantsTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MopUpRemnantsTool::Run(ThreeViewRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindBestShowers(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpRemnantsTool::FindBestShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
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

        TensorType::ElementList::const_iterator eIter = connectedElements.end();
        this->SelectBestElement(connectedElements, usedClusters, eIter);
        this->GetUsedClusters(connectedElements, usedClusters);

        if (connectedElements.end() == eIter)
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.emplace_back(eIter->GetClusterU());
        protoParticle.m_clusterList.emplace_back(eIter->GetClusterV());
        protoParticle.m_clusterList.emplace_back(eIter->GetClusterW());
        protoParticleVector.emplace_back(protoParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpRemnantsTool::GetUsedClusters(const TensorType::ElementList &elementList, ClusterSet &usedClusters) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        usedClusters.insert(eIter->GetClusterU());
        usedClusters.insert(eIter->GetClusterV());
        usedClusters.insert(eIter->GetClusterW());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpRemnantsTool::SelectBestElement(
    const TensorType::ElementList &elementList, const ClusterSet &usedClusters, TensorType::ElementList::const_iterator &bestIter) const
{
    float bestFigureOfMerit(0.f);

    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        const Cluster *const pClusterU = eIter->GetClusterU();
        const Cluster *const pClusterV = eIter->GetClusterV();
        const Cluster *const pClusterW = eIter->GetClusterW();

        if (usedClusters.count(pClusterU) || usedClusters.count(pClusterV) || usedClusters.count(pClusterW))
            continue;

        const float figureOfMerit(pClusterU->GetHadronicEnergy() + pClusterV->GetHadronicEnergy() + pClusterW->GetHadronicEnergy());

        if (figureOfMerit > bestFigureOfMerit)
        {
            bestFigureOfMerit = figureOfMerit;
            bestIter = eIter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MopUpRemnantsTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
