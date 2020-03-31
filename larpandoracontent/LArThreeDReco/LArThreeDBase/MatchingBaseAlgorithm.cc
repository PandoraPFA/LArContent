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

void MatchingBaseAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    if (!pInputClusterList)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    selectedClusterList = *pInputClusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::PrepareInputClusters(ClusterList &/*preparedClusterList*/)
{
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

void MatchingBaseAlgorithm::SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    this->SetPfoParticleId(pfoParameters);
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = protoParticle.m_clusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::SetPfoParticleId(PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    pfoParameters.m_particleId = E_MINUS; // Shower
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchingBaseAlgorithm::Run()
{
    try
    {
        this->SelectAllInputClusters();
        this->PrepareAllInputClusters();
        this->PerformMainLoop();
        this->ExamineOverlapContainer();
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
