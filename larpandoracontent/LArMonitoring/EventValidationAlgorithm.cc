/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_useTrueNeutrinosOnly(false),
    m_printAllToScreen(false),
    m_printMatchingToScreen(true),
    m_writeToTree(false),
    m_matchingMinSharedHits(5),
    m_matchingMinCompleteness(0.1f),
    m_matchingMinPurity(0.5f),
    m_fileIdentifier(0),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::~EventValidationAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::Run()
{
    ++m_eventNumber;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::MCContributionMap mcParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcParticleToHitsMap);

    if (!m_useTrueNeutrinosOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, mcParticleToHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, mcParticleToHitsMap);
    }

    PfoList finalStatePfos;
    for (const ParticleFlowObject *const pPfo : (pPfoList ? *pPfoList : PfoList()))
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }

    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, mcParticleToHitsMap, pfoToHitsMap);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {mcParticleToHitsMap}, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Print raw matching information to terminal
    if (m_printAllToScreen)
        this->PrintOutput(mcParticleToHitsMap, pfoToHitsMap, mcParticleToPfoHitSharingMap, false);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMCToPfoHitSharingMap(mcParticleToHitsMap, pfoToHitsMap, mcParticleToPfoHitSharingMap, interpretedMCToPfoHitSharingMap);

    if (m_printMatchingToScreen)
        this->PrintOutput(mcParticleToHitsMap, pfoToHitsMap, interpretedMCToPfoHitSharingMap, true);

//    if (m_writeToTree)
//        this->WriteOutput(mcParticleToHitsMap, pfoToHitsMap, interpretedMCToPfoHitSharingMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintOutput(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap, const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap,
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const bool printCorrectness) const
{
    std::cout << "---MATCHING-OUTPUT------------------------------------------------------------------------------" << std::endl;
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcParticleToHitsMap}, mcPrimaryVector);

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToHitsMap, primaryPfoVector);

    PfoToIdMap pfoToIdMap;
    unsigned int primaryIndex(0), pfoIndex(0), nCorrectNu(0), nCorrectTB(0), nCorrectCR(0), nTotalNu(0), nTotalTB(0), nTotalCR(0),
        nFakeNu(0), nSplitCR(0), nLostCR(0), nFakeCR(0);

    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
        (void) pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, pfoIndex++));

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const bool isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));
        const bool isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));

        if (isBeamNeutrinoFinalState) ++nTotalNu;
        if (isBeamParticle) ++nTotalTB;
        if (isCosmicRay) ++nTotalCR;

        const CaloHitList &mcPrimaryHitList(mcParticleToHitsMap.at(pMCPrimary));

        std::cout << std::endl << "--Primary " << primaryIndex++ // TODO nu, tb, cr index (nuance if applicable)
                  << ", Nu " << isBeamNeutrinoFinalState << ", TB " << isBeamParticle << ", CR " << isCosmicRay
                  << ", MCPDG " << pMCPrimary->GetParticleId()
                  << ", Energy " << pMCPrimary->GetEnergy()
                  << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                  << ", nMCHits " << mcPrimaryHitList.size()
                  << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList)
                  << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList)
                  << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        unsigned int nNuMatches(0), nCRMatches(0);

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcParticleToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(pfoToHitsMap.at(pfoToSharedHits.first));
            const bool isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));

            if (isRecoNeutrinoFinalState) ++nNuMatches;
            else ++nCRMatches;

            std::cout << "-MatchedPfo " << pfoToIdMap.at(pfoToSharedHits.first) // TODO nu, cr index
                << ", Nu " << isRecoNeutrinoFinalState << ", CR " << !isRecoNeutrinoFinalState
                << ", PDG " << pfoToSharedHits.first->GetParticleId()
                << ", nMatchedHits " << sharedHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                << ", nPfoHits " << pfoHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;
        }

        if (isBeamNeutrinoFinalState && (nNuMatches == 1) && (nCRMatches == 0)) ++nCorrectNu;
        else if (isBeamParticle && (nNuMatches == 1) && (nCRMatches == 0)) ++nCorrectTB;
        else if (isCosmicRay && (nNuMatches == 0) && (nCRMatches == 1)) ++nCorrectCR;
        else if (isCosmicRay && (nNuMatches > 0) && (nCRMatches == 0)) ++nFakeNu;
        else if (isCosmicRay && (nNuMatches + nCRMatches > 1)) ++nSplitCR;
        else if (isCosmicRay && (nNuMatches == 0) && (nCRMatches == 0)) ++nLostCR;
        else if (!isCosmicRay && (nCRMatches > 0)) ++nFakeCR;
    }

    if (printCorrectness)
    {
        std::cout << std::endl << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl
                  << "#Correct Nu primaries: " << nCorrectNu << ", Nu correct? " << (nTotalNu == nCorrectNu) << std::endl
                  << "#Correct TB particles: " << nCorrectTB << ", Fraction: " << (nTotalTB > 0 ? static_cast<float>(nCorrectTB) / static_cast<float>(nTotalTB) : 0.f) << std::endl
                  << "#Correct Cosmic Rays : " << nCorrectCR << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl
                  << "#CR as fake Nu: " << nFakeNu << ", #Split CRs: " << nSplitCR << ", #Lost CRs: " << nLostCR << ", #Fake CRs " << nFakeCR << std::endl;
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

//void EventValidationAlgorithm::WriteOutput(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
//    const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const
//{
//}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::InterpretMCToPfoHitSharingMap(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap,
    LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcParticleToHitsMap}, mcPrimaryVector);

    PfoSet usedPfos;
    while (this->GetStrongestPfoMatch(mcPrimaryVector, mcParticleToHitsMap, pfoToHitsMap, mcToPfoHitSharingMap, usedPfos, interpretedMCToPfoHitSharingMap)) {}
    this->GetRemainingPfoMatches(mcPrimaryVector, mcParticleToHitsMap, pfoToHitsMap, mcToPfoHitSharingMap, usedPfos, interpretedMCToPfoHitSharingMap);

    // Ensure all primaries have an entry, and sorting is as desired
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        LArMCParticleHelper::PfoToSharedHitsVector &pfoHitPairs(interpretedMCToPfoHitSharingMap[pMCPrimary]);
        std::sort(pfoHitPairs.begin(), pfoHitPairs.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool {
            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first)); });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::GetStrongestPfoMatch(const MCParticleVector &mcPrimaryVector, const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    const MCParticle *pBestMCParticle(nullptr);
    LArMCParticleHelper::PfoCaloHitListPair bestPfoHitPair(nullptr, CaloHitList());

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        // ATTN Used to have no constraints on input quality of mc primaries, then filter at this point
        if (interpretedMCToPfoHitSharingMap.count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            if (!this->IsGoodMatch(mcParticleToHitsMap.at(pMCPrimary), pfoToHitsMap.at(pfoToSharedHits.first), pfoToSharedHits.second))
                continue;

            if (pfoToSharedHits.second.size() > bestPfoHitPair.second.size())
            {
                pBestMCParticle = pMCPrimary;
                bestPfoHitPair = pfoToSharedHits;
            }
        }
    }

    if (!pBestMCParticle || !bestPfoHitPair.first)
        return false;

    interpretedMCToPfoHitSharingMap[pBestMCParticle].push_back(bestPfoHitPair);
    usedPfos.insert(bestPfoHitPair.first);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetRemainingPfoMatches(const MCParticleVector &mcPrimaryVector, const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap,
    const PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        // ATTN Used to have no constraints on input quality of mc primaries, then filter at this point
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            // ATTN Didn't used to include this here, when handling below threshold matches
            if (!this->IsGoodMatch(mcParticleToHitsMap.at(pMCPrimary), pfoToHitsMap.at(pfoToSharedHits.first), pfoToSharedHits.second))
                continue;

            const LArMCParticleHelper::MCParticleCaloHitListPair mcParticleToHits(pMCPrimary, pfoToSharedHits.second);
            LArMCParticleHelper::PfoToMCParticleHitSharingMap::iterator iter(pfoToMCParticleHitSharingMap.find(pfoToSharedHits.first));

            if (pfoToMCParticleHitSharingMap.end() == iter)
            {
                pfoToMCParticleHitSharingMap[pfoToSharedHits.first].push_back(mcParticleToHits);
            }
            else
            {
                if (1 != iter->second.size())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                LArMCParticleHelper::MCParticleCaloHitListPair &originalMCParticleToHits(iter->second.at(0));

                if (mcParticleToHits.second.size() > originalMCParticleToHits.second.size())
                    originalMCParticleToHits = mcParticleToHits;
            }
        }
    }

    for (const auto &mapEntry : pfoToMCParticleHitSharingMap)
    {
        const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleToHits(mapEntry.second.at(0));
        interpretedMCToPfoHitSharingMap[mcParticleToHits.first].push_back(LArMCParticleHelper::PfoCaloHitListPair(mapEntry.first, mcParticleToHits.second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsGoodMatch(const CaloHitList &trueHits, const CaloHitList &recoHits, const CaloHitList &sharedHits) const
{
    const float purity((recoHits.size() > 0) ? static_cast<float>(sharedHits.size()) / static_cast<float>(recoHits.size()) : 0.f);
    const float completeness((trueHits.size() > 0) ? static_cast<float>(sharedHits.size()) / static_cast<float>(trueHits.size()) : 0.f);

    return ((sharedHits.size() >= m_matchingMinSharedHits) && (purity >= m_matchingMinPurity) && (completeness >= m_matchingMinCompleteness));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
         "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinSharedHits", m_matchingMinSharedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinCompleteness", m_matchingMinCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPurity", m_matchingMinPurity));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "FileIdentifier", m_fileIdentifier));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
