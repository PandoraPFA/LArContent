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
    m_selectInputHits(true),
    m_minHitSharingFraction(0.9f),
    m_maxPhotonPropagation(2.5f),
    m_printAllToScreen(false),
    m_printMatchingToScreen(true),
    m_writeToTree(false),
    m_useSmallPrimaries(true),
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

    ValidationInfo validationInfo;
    this->FillValidationInfo(pMCParticleList, pCaloHitList, pPfoList, validationInfo);

    if (m_printAllToScreen)
        this->PrintOutput(validationInfo, false);

    if (m_printMatchingToScreen)
        this->PrintOutput(validationInfo, true);

    if (m_writeToTree)
        this->WriteOutput(validationInfo, true);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::PrimaryParameters parameters;

        parameters.m_selectInputHits = m_selectInputHits;
        parameters.m_minHitSharingFraction = m_minHitSharingFraction;
        parameters.m_maxPhotonPropagation = m_maxPhotonPropagation;
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }
    
    if (pPfoList)
    {
        PfoList finalStatePfos;
        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            if (LArPfoHelper::IsFinalState(pPfo))
                finalStatePfos.push_back(pPfo);
        }

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap);
        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching) const
{
    if (useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetAllMCParticleToHitsMap()}, mcPrimaryVector);

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    unsigned int primaryIndex(0), pfoIndex(0), nCorrectNu(0), nCorrectTB(0), nCorrectCR(0), nTotalNu(0), nTotalTB(0), nTotalCR(0),
        nFakeNu(0), nSplitCR(0), nLostCR(0), nFakeCR(0);

    PfoToIdMap pfoToIdMap;
    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, pfoIndex++));

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!hasMatch && !isTargetPrimary)
            continue;

        const bool isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const bool isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));
        const bool isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));

        if (isTargetPrimary && isBeamNeutrinoFinalState) ++nTotalNu;
        if (isTargetPrimary && isBeamParticle) ++nTotalTB;
        if (isTargetPrimary && isCosmicRay) ++nTotalCR;

        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        std::cout << std::endl << (!isTargetPrimary ? "(Non target) " : "")
                  << "Primary " << primaryIndex++
                  << ", Nu " << isBeamNeutrinoFinalState << ", TB " << isBeamParticle << ", CR " << isCosmicRay
                  << ", MCPDG " << pMCPrimary->GetParticleId()
                  << ", Energy " << pMCPrimary->GetEnergy()
                  << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                  << ", nMCHits " << mcPrimaryHitList.size()
                  << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList)
                  << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList)
                  << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        unsigned int nNuMatches(0), nCRMatches(0);

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            if (isGoodMatch && isRecoNeutrinoFinalState) ++nNuMatches;
            else if (isGoodMatch) ++nCRMatches;

            std::cout << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                      << "MatchedPfo " << pfoToIdMap.at(pfoToSharedHits.first)
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

        if (isTargetPrimary && isBeamNeutrinoFinalState && (nNuMatches == 1) && (nCRMatches == 0)) ++nCorrectNu;
        else if (isTargetPrimary && isBeamParticle && (nNuMatches == 1) && (nCRMatches == 0)) ++nCorrectTB;
        else if (isTargetPrimary && isCosmicRay && (nNuMatches == 0) && (nCRMatches == 1)) ++nCorrectCR;
        else if (isTargetPrimary && isCosmicRay && (nNuMatches > 0) && (nCRMatches == 0)) ++nFakeNu;
        else if (isTargetPrimary && isCosmicRay && (nNuMatches + nCRMatches > 1)) ++nSplitCR;
        else if (isTargetPrimary && isCosmicRay && (nNuMatches == 0) && (nCRMatches == 0)) ++nLostCR;
        else if (isTargetPrimary && !isCosmicRay && (nCRMatches > 0)) ++nFakeCR;
    }

    if (useInterpretedMatching)
    {
        std::cout << std::endl << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalNu > 0) std::cout << "#Correct target Nu primaries: " << nCorrectNu << "/" << nTotalNu << ", Nu correct? " << (nTotalNu == nCorrectNu) << std::endl;
        if (nTotalTB > 0) std::cout << "#Correct target TB particles: " << nCorrectTB << "/" << nTotalTB << ", Fraction: " << (nTotalTB > 0 ? static_cast<float>(nCorrectTB) / static_cast<float>(nTotalTB) : 0.f) << std::endl;
        if (nTotalCR > 0) std::cout << "#Correct target Cosmic Rays : " << nCorrectCR << "/" << nTotalCR << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if ((nTotalCR > 0) || (nFakeCR > 0)) std::cout << "#CR as fake Nu: " << nFakeNu << ", #Split CRs: " << nSplitCR << ", #Lost CRs: " << nLostCR << ", #Fake CRs " << nFakeCR << std::endl;;
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteOutput(const ValidationInfo &/*validationInfo*/, const bool /*useInterpretedMatching*/) const
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::InterpretMatching(const ValidationInfo &validationInfo, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetAllMCParticleToHitsMap()}, mcPrimaryVector);

    PfoSet usedPfos;
    while (this->GetStrongestPfoMatch(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap)) {}
    this->GetRemainingPfoMatches(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap);

    // Ensure all primaries have an entry, and sorting is as desired
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        LArMCParticleHelper::PfoToSharedHitsVector &pfoHitPairs(interpretedMCToPfoHitSharingMap[pMCPrimary]);
        std::sort(pfoHitPairs.begin(), pfoHitPairs.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool {
            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first)); });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::GetStrongestPfoMatch(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    const MCParticle *pBestMCParticle(nullptr);
    LArMCParticleHelper::PfoCaloHitListPair bestPfoHitPair(nullptr, CaloHitList());

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (interpretedMCToPfoHitSharingMap.count(pMCPrimary))
            continue;

        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            if (!this->IsGoodMatch(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary), validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first), pfoToSharedHits.second))
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

void EventValidationAlgorithm::GetRemainingPfoMatches(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    const PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseSmallPrimaries", m_useSmallPrimaries));

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
