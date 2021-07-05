/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the cheating neutrino creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoHierarchyAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoHierarchyAlgorithm::CheatingNeutrinoHierarchyAlgorithm() : m_collapseToPrimaryMCParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoHierarchyAlgorithm::Run()
{
    const ParticleFlowObject *pNeutrinoPfo(nullptr);
    this->GetNeutrinoPfo(pNeutrinoPfo);

    MCParticleVector mcNeutrinoVector;
    this->GetMCNeutrinoVector(mcNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        MCParticleToPfoMap mcParticleToPfoMap;
        this->GetMCParticleToDaughterPfoMap(mcParticleToPfoMap);
        this->CreatePfoHierarchy(pMCNeutrino, pNeutrinoPfo, mcParticleToPfoMap);
    }

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoHierarchyAlgorithm::GetNeutrinoPfo(const ParticleFlowObject *&pNeutrinoPfo) const
{
    const PfoList *pPfoList(nullptr);

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
	  std::cout << "CheatingNeutrinoHierarchyAlgorithm: unable to find pfo list " << m_neutrinoPfoListName << std::endl;

        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    // ATTN No way of linking pfo neutrino to MC neutrino so enforces that only one pfo, of neutrino-type, be in the specified input list
    pNeutrinoPfo = ((1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr);

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoHierarchyAlgorithm::GetMCNeutrinoVector(MCParticleVector &mcNeutrinoVector) const
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    MCParticleVector allMCNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, allMCNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : allMCNeutrinoVector)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCNeutrino) && (MC_3D == pMCNeutrino->GetMCParticleType()))
            mcNeutrinoVector.push_back(pMCNeutrino);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoHierarchyAlgorithm::GetMCParticleToDaughterPfoMap(MCParticleToPfoMap &mcParticleToPfoMap) const
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;

    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }

    for (const std::string &daughterPfoListName : m_daughterPfoListNames)
    {
        const PfoList *pDaughterPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, daughterPfoListName, pDaughterPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CheatingNeutrinoHierarchyAlgorithm: pfo list " << daughterPfoListName << " unavailable." << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pDaughterPfo : *pDaughterPfoList)
        {
            const MCParticle *pMCParticle(LArMCParticleHelper::GetMainMCParticle(pDaughterPfo));

            if (m_collapseToPrimaryMCParticles)
            {
                LArMCParticleHelper::MCRelationMap::const_iterator primaryIter = mcPrimaryMap.find(pMCParticle);

                if (mcPrimaryMap.end() == primaryIter)
                    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                pMCParticle = primaryIter->second;
            }

            if (!mcParticleToPfoMap.insert(MCParticleToPfoMap::value_type(pMCParticle, pDaughterPfo)).second)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoHierarchyAlgorithm::CreatePfoHierarchy(const MCParticle *const pParentMCParticle,
    const ParticleFlowObject *const pParentPfo, const MCParticleToPfoMap &mcParticleToPfoMap) const
{
    for (const MCParticle *const pDaughterMCParticle : pParentMCParticle->GetDaughterList())
    {
        MCParticleToPfoMap::const_iterator mapIter = mcParticleToPfoMap.find(pDaughterMCParticle);

        if (mcParticleToPfoMap.end() == mapIter)
        {
            this->CreatePfoHierarchy(pDaughterMCParticle, pParentPfo, mcParticleToPfoMap);
            continue;
        }

        const ParticleFlowObject *const pDaughterPfo(mapIter->second);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));
        this->CreatePfoHierarchy(pDaughterMCParticle, pDaughterPfo, mcParticleToPfoMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterPfoListNames", m_daughterPfoListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
