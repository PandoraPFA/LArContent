/**
 *  @file   LArContent/src/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArMonitoringHelper.h"

#include "LArMonitoring/EventValidationAlgorithm.h"

#include "LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_extractNeutrinoDaughters(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::~EventValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::Run()
{
#ifdef MONITORING
    // Input collections
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = NULL;
    PfoList pfoList((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) ? PfoList(*pPfoList) : PfoList());

    if (m_extractNeutrinoDaughters)
        LArMonitoringHelper::ExtractNeutrinoDaughters(pfoList);

    // Extract monitoring information
    MCParticleVector mcNeutrinoList;                                // true neutrinos
    LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoList);

    MCParticleVector mcPrimaryList;                                 // primary mc particles
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryList);

    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;            // [mc particles -> primary mc particle]
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    LArMonitoringHelper::CaloHitToMCMap hitToPrimaryMCMap;          // [hit -> primary mc particle]
    LArMonitoringHelper::MCContributionMap mcToTrueHitListMap;      // [primary mc particle -> true hit list]
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

    LArMonitoringHelper::CaloHitToPfoMap hitToPfoMap;               // [hit -> pfo]
    LArMonitoringHelper::PfoContributionMap pfoToHitListMap;        // [pfo -> reco hit list]
    LArMonitoringHelper::GetPfoToCaloHitMatches(pCaloHitList, pfoList, hitToPfoMap, pfoToHitListMap);

    LArMonitoringHelper::MCToPfoMap mcToBestPfoMap;                 // [mc particle -> best matched pfo]
    LArMonitoringHelper::MCContributionMap mcToBestPfoHitsMap;      // [mc particle -> list of hits included in best pfo]
    LArMonitoringHelper::MCToPfoMatchingMap mcToFullPfoMatchingMap; // [mc particle -> all matched pfos (and matched hits)]
    LArMonitoringHelper::GetMCParticleToPfoMatches(pCaloHitList, pfoList, hitToPrimaryMCMap, mcToBestPfoMap, mcToBestPfoHitsMap, mcToFullPfoMatchingMap);

    // Use monitoring information
    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);

        if (!pLArMCNeutrino)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        std::cout << "MCNeutrino " << pLArMCNeutrino << ", PDG " << pLArMCNeutrino->GetParticleId() << ", Nuance " << pLArMCNeutrino->GetNuanceCode() << std::endl;
    }

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        std::cout << "-Primary " << pMCPrimary << ", PDG " << pMCPrimary->GetParticleId() << std::endl;

        LArMonitoringHelper::MCContributionMap::const_iterator trueHitsIter = mcToTrueHitListMap.find(pMCPrimary);

        if (mcToTrueHitListMap.end() == trueHitsIter)
        {
            std::cout << "--No true hits " << std::endl;
        }
        else 
        {
            const CaloHitList &caloHitList(trueHitsIter->second);
            std::cout << "--True hits " << caloHitList.size()
                      << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList)
                      << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList)
                      << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList) << ") " << std::endl;
        }

        LArMonitoringHelper::MCToPfoMatchingMap::const_iterator matchedPfoIter = mcToFullPfoMatchingMap.find(pMCPrimary);

        if ((mcToFullPfoMatchingMap.end() == matchedPfoIter) || matchedPfoIter->second.empty())
        {
            std::cout << "--No matched pfo " << std::endl;
        }
        else
        {
            for (const LArMonitoringHelper::PfoContributionMap::value_type contribution : matchedPfoIter->second)
            {
                const ParticleFlowObject *const pMatchedPfo(contribution.first);
                const CaloHitList &matchedCaloHitList(contribution.second);

                LArMonitoringHelper::PfoContributionMap::const_iterator pfoHitsIter = pfoToHitListMap.find(pMatchedPfo);

                if (pfoToHitListMap.end() == pfoHitsIter)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const CaloHitList &pfoCaloHitList(pfoHitsIter->second);

                std::cout << "--MatchedPfo " << pMatchedPfo << ", PDG " << pMatchedPfo->GetParticleId()
                          << ", pfo hits " << pfoCaloHitList.size()
                          << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoCaloHitList)
                          << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoCaloHitList)
                          << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoCaloHitList) << ") "
                          << ", matched hits " << matchedCaloHitList.size()
                          << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, matchedCaloHitList)
                          << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, matchedCaloHitList)
                          << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, matchedCaloHitList) << ") " << std::endl;
            }
        }
    }

    // Write monitoring tree
    IntVector nMatchedHitsWVector;
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatchedHitsW", &nMatchedHitsWVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
#endif
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtractNeutrinoDaughters", m_extractNeutrinoDaughters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
