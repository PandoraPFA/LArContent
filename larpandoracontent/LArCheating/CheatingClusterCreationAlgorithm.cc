/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterCreationAlgorithm::CheatingClusterCreationAlgorithm() :
    m_collapseToPrimaryMCParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCreationAlgorithm::Run()
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;

    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    MCParticleToHitListMap mcParticleToHitListMap;
    this->GetMCParticleToHitListMap(pCaloHitList, mcPrimaryMap, mcParticleToHitListMap);

    this->CreateClusters(mcParticleToHitListMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingClusterCreationAlgorithm::GetMCParticleToHitListMap(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap,
    MCParticleToHitListMap &mcParticleToHitListMap) const
{
    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            if (!PandoraContentApi::IsAvailable(*this, *iter))
                continue;

            this->SimpleMCParticleCollection(*iter, mcPrimaryMap, mcParticleToHitListMap);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingClusterCreationAlgorithm::SimpleMCParticleCollection(const CaloHit *const pCaloHit, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap,
    MCParticleToHitListMap &mcParticleToHitListMap) const
{
    const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

    if (!this->SelectMCParticlesForClustering(pMCParticle))
        return;

    if (m_collapseToPrimaryMCParticles)
    {
        LArMCParticleHelper::MCRelationMap::const_iterator primaryIter = mcPrimaryMap.find(pMCParticle);

        if (mcPrimaryMap.end() == primaryIter)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        pMCParticle = primaryIter->second;
    }

    mcParticleToHitListMap[pMCParticle].push_back(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingClusterCreationAlgorithm::SelectMCParticlesForClustering(const MCParticle *const pMCParticle) const
{
    if (m_particleIdList.empty())
        return true;

    for (IntVector::const_iterator iter = m_particleIdList.begin(), iterEnd = m_particleIdList.end(); iter != iterEnd; ++iter)
    {
        if (pMCParticle->GetParticleId() == *iter)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingClusterCreationAlgorithm::CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const
{
    MCParticleVector mcParticleVector;
    for (const auto &mapEntry : mcParticleToHitListMap) mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const CaloHitList &caloHitList(mcParticleToHitListMap.at(pMCParticle));

        if (caloHitList.empty())
            continue;

        const Cluster *pCluster = NULL;
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

        PandoraContentApi::Cluster::Metadata metadata;

        switch (pMCParticle->GetParticleId())
        {
        case PHOTON:
        case E_PLUS:
        case E_MINUS:
        case MU_PLUS:
        case MU_MINUS:
            metadata.m_particleId = pMCParticle->GetParticleId();
            break;
        default:
            break;
        }

        if (metadata.m_particleId.IsInitialized())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    if (m_collapseToPrimaryMCParticles)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "MCParticleListName", m_mcParticleListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ParticleIdList", m_particleIdList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
