/**
 *  @file   larpandoracontent/LArCheating/CheatingCCLeptonClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingCCLeptonClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingCCLeptonClusterCreationAlgorithm::CheatingCCLeptonClusterCreationAlgorithm() : m_collapseToPrimaryMCParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCLeptonClusterCreationAlgorithm::Run()
{
    MCParticleToHitListMap mcParticleToHitListMap;
    this->GetMCParticleToHitListMap(mcParticleToHitListMap);
    this->CreateClusters(mcParticleToHitListMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::GetMCParticleToHitListMap(MCParticleToHitListMap &mcParticleToHitListMap) const
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;

    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            this->SimpleMCParticleCollection(pCaloHit, mcPrimaryMap, mcParticleToHitListMap);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::SimpleMCParticleCollection(const CaloHit *const pCaloHit,
    const LArMCParticleHelper::MCRelationMap &mcPrimaryMap, MCParticleToHitListMap &mcParticleToHitListMap) const
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

bool CheatingCCLeptonClusterCreationAlgorithm::SelectMCParticlesForClustering(const MCParticle *const pMCParticle) const
{
    // input def of is this a leading CC lepton?
    // this might fail?
    std::cout << "ISOBEL YOU STILL NEED TO IMPLEMENT WHETHER THIS IS A CC OR NC INTERACTION" << std::endl;

    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));
    const int pdg(pParentMCParticle->GetParticleId());
    const bool isLeadingLepton((E_MINUS == std::abs(pdg)) || (MU_MINUS == std::abs(pdg)));

    return isLeadingLepton;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const
{
    MCParticleVector mcParticleVector;
    for (const auto &mapEntry : mcParticleToHitListMap)
        mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const CaloHitList &caloHitList(mcParticleToHitListMap.at(pMCParticle));

        if (caloHitList.empty())
            continue;

        const Cluster *pCluster(nullptr);
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

StatusCode CheatingCCLeptonClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    if (m_collapseToPrimaryMCParticles)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
