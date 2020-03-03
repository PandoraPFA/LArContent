/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.cc
 *
 *  @brief  Implementation of the three dimension algorithm base class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MatchingBaseAlgorithm::MatchingBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MatchingBaseAlgorithm::~MatchingBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::CreateThreeDParticles(const ProtoParticleVector &protoParticleVector)
{
    bool particlesMade(false);
    const PfoList *pPfoList(nullptr); std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const ProtoParticle &protoParticle : protoParticleVector)
    {
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        this->SetPfoParameters(protoParticle, pfoParameters);

        const ParticleFlowObject *pPfo(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
        particlesMade = true;
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::MakeClusterMerges(const ClusterMergeMap &clusterMergeMap)
{
    ClusterSet deletedClusters;

    ClusterList parentClusters;
    for (const auto &mapEntry : clusterMergeMap) parentClusters.push_back(mapEntry.first);
    parentClusters.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : parentClusters)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pParentCluster));
        const std::string &clusterListName(this->GetClusterListName(hitType));

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterList daughterClusters(clusterMergeMap.at(pParentCluster));
        daughterClusters.sort(LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : daughterClusters)
        {
            if (deletedClusters.count(pParentCluster) || deletedClusters.count(pDaughterCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            this->UpdateUponDeletion(pDaughterCluster);
            this->UpdateUponDeletion(pParentCluster);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster, clusterListName, clusterListName));

            this->UpdateForNewCluster(pParentCluster);
            deletedClusters.insert(pDaughterCluster);
        }
    }

    return !(deletedClusters.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    // TODO Correct these placeholder parameters
    pfoParameters.m_particleId = E_MINUS; // Shower
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListU.begin(), protoParticle.m_clusterListU.end());
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListV.begin(), protoParticle.m_clusterListV.end());
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListW.begin(), protoParticle.m_clusterListW.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::SelectAllInputClusters()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::PreparationStep()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::TidyUp()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchingBaseAlgorithm::Run()
{
    try
    {
        this->SelectAllInputClusters();
        this->PreparationStep();
        this->PerformMainLoop();
        this->ExamineTensor();
        this->TidyUp();
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();

        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
