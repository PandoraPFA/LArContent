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
#include "LArHelpers/LArPfoHelper.h"

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
    PfoIdMap pfoIdMap;                                              // pfo -> unique identifier
    this->GetPfoIdMap(pfoList, pfoIdMap);

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
    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);

        if (!pLArMCNeutrino)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        std::cout << "MCNeutrino, PDG " << pLArMCNeutrino->GetParticleId() << ", Nuance " << pLArMCNeutrino->GetNuanceCode() << std::endl;
    }

    SimpleMCPrimaryList simpleMCPrimaryList;

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        SimpleMCPrimary simpleMCPrimary;
        simpleMCPrimary.m_pPandoraAddress = pMCPrimary;
        simpleMCPrimary.m_pdgCode = pMCPrimary->GetParticleId();

        LArMonitoringHelper::MCContributionMap::const_iterator trueHitsIter = mcToTrueHitListMap.find(pMCPrimary);

        if (mcToTrueHitListMap.end() != trueHitsIter)
        {
            const CaloHitList &caloHitList(trueHitsIter->second);
            simpleMCPrimary.m_nMCHitsTotal = caloHitList.size();
            simpleMCPrimary.m_nMCHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList);
            simpleMCPrimary.m_nMCHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList);
            simpleMCPrimary.m_nMCHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList);
        }

        LArMonitoringHelper::MCToPfoMatchingMap::const_iterator matchedPfoIter = mcToFullPfoMatchingMap.find(pMCPrimary);

        if (mcToFullPfoMatchingMap.end() != matchedPfoIter)
            simpleMCPrimary.m_nMatchedPfos = matchedPfoIter->second.size();

        simpleMCPrimaryList.push_back(simpleMCPrimary);
    }

    int primaryId(0);
    std::sort(simpleMCPrimaryList.begin(), simpleMCPrimaryList.end(), EventValidationAlgorithm::SortSimpleMCPrimaries);

    for (const SimpleMCPrimary &simpleMCPrimary : simpleMCPrimaryList)
    {
        std::cout << std::endl << "Primary " << primaryId++ << ", PDG " << simpleMCPrimary.m_pdgCode
                  << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << ")"
                  << std::endl;

        SimpleMatchedPfoList simpleMatchedPfoList;
        LArMonitoringHelper::MCToPfoMatchingMap::const_iterator matchedPfoIter = mcToFullPfoMatchingMap.find(simpleMCPrimary.m_pPandoraAddress);

        if (mcToFullPfoMatchingMap.end() != matchedPfoIter)
        {
            for (const LArMonitoringHelper::PfoContributionMap::value_type contribution : matchedPfoIter->second)
            {
                const ParticleFlowObject *const pMatchedPfo(contribution.first);
                const CaloHitList &matchedCaloHitList(contribution.second);

                SimpleMatchedPfo simpleMatchedPfo;
                simpleMatchedPfo.m_pPandoraAddress = pMatchedPfo;
                simpleMatchedPfo.m_id = pfoIdMap.at(pMatchedPfo);
                simpleMatchedPfo.m_pdgCode = pMatchedPfo->GetParticleId();

                simpleMatchedPfo.m_nMatchedHitsTotal = matchedCaloHitList.size();
                simpleMatchedPfo.m_nMatchedHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, matchedCaloHitList);
                simpleMatchedPfo.m_nMatchedHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, matchedCaloHitList);
                simpleMatchedPfo.m_nMatchedHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, matchedCaloHitList);

                LArMonitoringHelper::PfoContributionMap::const_iterator pfoHitsIter = pfoToHitListMap.find(pMatchedPfo);

                if (pfoToHitListMap.end() == pfoHitsIter)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const CaloHitList &pfoCaloHitList(pfoHitsIter->second);

                simpleMatchedPfo.m_nPfoHitsTotal = pfoCaloHitList.size();
                simpleMatchedPfo.m_nPfoHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoCaloHitList);
                simpleMatchedPfo.m_nPfoHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoCaloHitList);
                simpleMatchedPfo.m_nPfoHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoCaloHitList);

                simpleMatchedPfoList.push_back(simpleMatchedPfo);
            }
        }

        std::sort(simpleMatchedPfoList.begin(), simpleMatchedPfoList.end(), EventValidationAlgorithm::SortSimpleMatchedPfos);

        for (const SimpleMatchedPfo simpleMatchedPfo : simpleMatchedPfoList)
        {
            std::cout << "-MatchedPfo " << simpleMatchedPfo.m_id << ", PDG " << simpleMatchedPfo.m_pdgCode
                      << ", nMatchedHits " << simpleMatchedPfo.m_nMatchedHitsTotal << " (" << simpleMatchedPfo.m_nMatchedHitsU << ", " << simpleMatchedPfo.m_nMatchedHitsV << ", " << simpleMatchedPfo.m_nMatchedHitsW << ")"
                      << ", nPfoHits " << simpleMatchedPfo.m_nPfoHitsTotal << " (" << simpleMatchedPfo.m_nPfoHitsU << ", " << simpleMatchedPfo.m_nPfoHitsV << ", " << simpleMatchedPfo.m_nPfoHitsW << ")"
                      << std::endl;
        }
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    // Write monitoring tree
    IntVector nMatchedHitsWVector;
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatchedHitsW", &nMatchedHitsWVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
#endif
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetPfoIdMap(const pandora::PfoList &pfoList, PfoIdMap &pfoIdMap) const
{
    PfoVector pfoVector(pfoList.begin(), pfoList.end());
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    int id(0);

    for (const ParticleFlowObject *const pPfo : pfoVector)
        (void) pfoIdMap.insert(PfoIdMap::value_type(pPfo, id++));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::SortSimpleMCPrimaries(const SimpleMCPrimary &lhs, const SimpleMCPrimary &rhs)
{
    return (lhs.m_nMCHitsTotal > rhs.m_nMCHitsTotal);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::SortSimpleMatchedPfos(const SimpleMatchedPfo &lhs, const SimpleMatchedPfo &rhs)
{
    return (lhs.m_nMatchedHitsTotal > rhs.m_nMatchedHitsTotal);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::SimpleMCPrimary::SimpleMCPrimary() :
    m_pdgCode(0),
    m_nMCHitsTotal(0),  
    m_nMCHitsU(0),
    m_nMCHitsV(0),
    m_nMCHitsW(0),
    m_nMatchedPfos(0),
    m_pPandoraAddress(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::SimpleMatchedPfo::SimpleMatchedPfo() :
    m_id(-1),
    m_pdgCode(0), 
    m_nPfoHitsTotal(0),
    m_nPfoHitsU(0),
    m_nPfoHitsV(0),
    m_nPfoHitsW(0),
    m_nMatchedHitsTotal(0),
    m_nMatchedHitsU(0),
    m_nMatchedHitsV(0),
    m_nMatchedHitsW(0),
    m_pPandoraAddress(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
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
