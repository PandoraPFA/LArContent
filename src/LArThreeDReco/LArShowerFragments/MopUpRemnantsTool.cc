/**
 *  @file   LArContent/src/LArThreeDReco/ShowerFragments/MopUpRemnantsTool.cc
 *
 *  @brief  Implementation of the mop-up remnants tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArShowerFragments/MopUpRemnantsTool.h"

using namespace pandora;

namespace lar_content
{

MopUpRemnantsTool::MopUpRemnantsTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MopUpRemnantsTool::Run(ThreeDRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindBestShowers(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpRemnantsTool::FindBestShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        TensorType::ElementList connectedElements;
        overlapTensor.GetConnectedElements(iterU->first, true, connectedElements);

        TensorType::ElementList::const_iterator eIter = connectedElements.end();
        this->SelectBestElement(connectedElements, usedClusters, eIter);
        this->GetClusters(connectedElements, usedClusters);

        if (connectedElements.end() == eIter)
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(eIter->GetClusterU());
        protoParticle.m_clusterListV.insert(eIter->GetClusterV());
        protoParticle.m_clusterListW.insert(eIter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpRemnantsTool::GetClusters(const TensorType::ElementList &elementList, ClusterList &clusterList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        clusterList.insert(eIter->GetClusterU());
        clusterList.insert(eIter->GetClusterV());
        clusterList.insert(eIter->GetClusterW());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpRemnantsTool::SelectBestElement(const TensorType::ElementList &elementList, const ClusterList &usedClusters,
    TensorType::ElementList::const_iterator &bestIter) const
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
