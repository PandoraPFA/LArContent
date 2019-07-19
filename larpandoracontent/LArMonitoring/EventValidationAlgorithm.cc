/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_useTrueNeutrinosOnly(false),
    m_testBeamMode(false),
    m_testBeamHierarchyMode(false),
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
        this->PrintAllMatches(validationInfo);

    if (m_printMatchingToScreen)
        this->PrintInterpretedMatches(validationInfo);

    if (m_writeToTree)
        this->WriteInterpretedMatches(validationInfo);

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
        if (!m_useTrueNeutrinosOnly && !m_testBeamHierarchyMode) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);
        if (m_testBeamHierarchyMode) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsLeadingBeamParticle, targetMCParticleToHitsMap, true);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly && !m_testBeamHierarchyMode) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, allMCParticleToHitsMap);
        if (m_testBeamHierarchyMode) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsLeadingBeamParticle, allMCParticleToHitsMap, true);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {
        PfoList allConnectedPfos;
        LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

        PfoList finalStatePfos;
        for (const ParticleFlowObject *const pPfo : allConnectedPfos)
        {
            // ATTN: Is test beam only set for parent pfo, therefor add parent and daughters for that particle
            if (m_testBeamHierarchyMode)
            {
                if (LArPfoHelper::IsTestBeam(pPfo))
                {
                    finalStatePfos.push_back(pPfo);
                    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
                        finalStatePfos.push_back(pDaughterPfo);
                }
                else if (pPfo->GetParentPfoList().empty())
                {
                    finalStatePfos.push_back(pPfo);
                }
            }
            else
            {
                if ((!m_testBeamMode && LArPfoHelper::IsFinalState(pPfo)) || (m_testBeamMode && pPfo->GetParentPfoList().empty()))
                    finalStatePfos.push_back(pPfo);
            }
        }

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;

        if (m_testBeamHierarchyMode)
        {
            LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, true);
        }
        else
        {
            LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap);
        }

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

void EventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    // Neutrino Validation Bookkeeping
    int nNeutrinoPrimaries(0);
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary)) ++nNeutrinoPrimaries;

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    int pfoIndex(0), neutrinoPfoIndex(0), testBeamPfoIndex(0);
    PfoToIdMap pfoToIdMap, neutrinoPfoToIdMap, testBeamPfoToIdMap;

    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));

        const Pfo *const pRecoNeutrino(LArPfoHelper::IsNeutrinoFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentNeutrino(pPrimaryPfo) : nullptr);
        const Pfo *const pRecoTestBeam(LArPfoHelper::IsTestBeamFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentPfo(pPrimaryPfo) : nullptr);

        if (pRecoNeutrino && !neutrinoPfoToIdMap.count(pRecoNeutrino))
            neutrinoPfoToIdMap.insert(PfoToIdMap::value_type(pRecoNeutrino, ++neutrinoPfoIndex));

        if (pRecoTestBeam && !testBeamPfoToIdMap.count(pRecoTestBeam))
            testBeamPfoToIdMap.insert(PfoToIdMap::value_type(pRecoTestBeam, ++testBeamPfoIndex));
    }

    LArMCParticleHelper::MCParticleIntMap triggeredToLeading, triggeredToLeadingCounter;

    // ProtoDUNE Hierarchy Validation Bookkeeping
    if (m_testBeamHierarchyMode)
    {
        // ATTN: At this stage the mcPrimaryVector is ordered from neutrinos, beam and then cosmics.  Here we extract the beam and reorder to ensure
        // the order follows primary parent beam 1, daughter 1 of beam 1, daughter 2 of beam 1, ..., primary parent beam 2, daughter 1 of beam 2, etc...
        // as expected by downstream logic
        MCParticleVector temporaryMCParticleVector(mcPrimaryVector), triggeredBeamParticles;
        LArMCParticleHelper::MCRelationMap leadingToTriggeredMap;
        mcPrimaryVector.clear();

        for (const MCParticle *const pMCPrimary : temporaryMCParticleVector)
        {
            if (LArMCParticleHelper::IsLeadingBeamParticle(pMCPrimary))
            {
                const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCPrimary));
                leadingToTriggeredMap.insert(LArMCParticleHelper::MCRelationMap::value_type(pMCPrimary, pParentMCParticle));

                if (std::find(triggeredBeamParticles.begin(), triggeredBeamParticles.end(), pParentMCParticle) == triggeredBeamParticles.end())
                    triggeredBeamParticles.push_back(pParentMCParticle);
            }
            else
            {
                mcPrimaryVector.push_back(pMCPrimary);
            }
        }

        for (const MCParticle *const pMCParent : triggeredBeamParticles)
        {
            // Parent appears first
            mcPrimaryVector.push_back(pMCParent);
            triggeredToLeading.insert(LArMCParticleHelper::MCParticleIntMap::value_type(pMCParent, 1));
            triggeredToLeadingCounter.insert(LArMCParticleHelper::MCParticleIntMap::value_type(pMCParent, 0));

            for (const auto iter : leadingToTriggeredMap)
            {
                // Followed by daughters, veto parent <-> parent matche
                if (iter.second == pMCParent && iter.first != pMCParent)
                {
                    mcPrimaryVector.push_back(iter.first);
                    triggeredToLeading.at(pMCParent)++;
                }
            }
        }
    }

    PfoSet recoNeutrinos, recoTestBeamHierarchies;
    MCParticleList associatedMCPrimaries;

    int nCorrectNu(0), nTotalNu(0), nCorrectTB(0), nTotalTB(0), nCorrectTBHierarchy(0), nTotalTBHierarchy(0), nCorrectCR(0), nTotalCR(0), nFakeNu(0), nFakeCR(0), nSplitNu(0), nSplitCR(0), nLost(0);
    int mcPrimaryIndex(0), nTargetMatches(0), nTargetNuMatches(0), nTargetCRMatches(0), nTargetGoodNuMatches(0), nTargetNuSplits(0), nTargetNuLosses(0);
    IntVector mcPrimaryId, mcPrimaryPdg, mcPrimaryTier, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedNuPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoTier, bestMatchPfoIsRecoNu, bestMatchPfoRecoNuId, bestMatchPfoIsTestBeam, bestMatchPfoRecoTBId;
    IntVector bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;

    std::stringstream targetSS;
    std::string name(m_testBeamHierarchyMode ? "TB" : "Nu");

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!m_testBeamHierarchyMode)
        {
            if (!isTargetPrimary && !hasMatch)
                continue;
        }
        else
        {
            if (!isTargetPrimary)
                continue;

            // Parent in hierarchy needed even if no match
            const bool hasVisibleTargets(!triggeredToLeading.empty() ? triggeredToLeading.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)) != 1 : false);

            if (!hasMatch && !hasVisibleTargets)
                continue;
        }

        associatedMCPrimaries.push_back(pMCPrimary);
        const int nTargetPrimaries(associatedMCPrimaries.size());
        const bool isLastNeutrinoPrimary(++mcPrimaryIndex == nNeutrinoPrimaries);
        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        const int isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const int isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));

        // Leading beam particle is the primary beam particle or a daughter of that particle
        const int isLeadingBeamParticle(LArMCParticleHelper::IsLeadingBeamParticle(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));

        // Tier (0) : Primary, (1) : Daughter, (2) : Granddaughter etc...  Note that the tier increases for both visible and invisible particles
        const int mcHierarchyTier(LArMCParticleHelper::GetHierarchyTier(pMCPrimary));

        // Identify the number of matched leading particles and flag whether last particle in hierarchy is being considered
        bool isLastTestBeamLeading(false);
        if (isLeadingBeamParticle && m_testBeamHierarchyMode)
        {
            triggeredToLeadingCounter.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary))++;
            const int nHierarchyLeading(triggeredToLeadingCounter.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
            isLastTestBeamLeading = (nHierarchyLeading == triggeredToLeading.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        }

#ifdef MONITORING
        const CartesianVector &targetVertex(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetVertex());
        const float targetVertexX(targetVertex.GetX()), targetVertexY(targetVertex.GetY()), targetVertexZ(targetVertex.GetZ());
#endif

        if (m_testBeamHierarchyMode) for (int tier = 0; tier < mcHierarchyTier; tier++) targetSS << " -> ";

        targetSS << (!isTargetPrimary ? "(Non target) " : "")
                 << "PrimaryId " << mcPrimaryIndex
                 << ", Nu " << isBeamNeutrinoFinalState
                 << ", TB " << isBeamParticle;
        if (m_testBeamHierarchyMode) targetSS << ", TB Hierarchy " << isLeadingBeamParticle;
        targetSS << ", CR " << isCosmicRay
                 << ", MCPDG " << pMCPrimary->GetParticleId();
        if (m_testBeamHierarchyMode) targetSS << ", Tier " << mcHierarchyTier;
        targetSS << ", Energy " << pMCPrimary->GetEnergy()
                 << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                 << ", nMCHits " << mcPrimaryHitList.size()
                 << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        mcPrimaryId.push_back(mcPrimaryIndex);
        mcPrimaryPdg.push_back(pMCPrimary->GetParticleId());
        mcPrimaryTier.push_back(mcHierarchyTier);
        mcPrimaryE.push_back(pMCPrimary->GetEnergy());
        mcPrimaryPX.push_back(pMCPrimary->GetMomentum().GetX());
        mcPrimaryPY.push_back(pMCPrimary->GetMomentum().GetY());
        mcPrimaryPZ.push_back(pMCPrimary->GetMomentum().GetZ());
        mcPrimaryVtxX.push_back(pMCPrimary->GetVertex().GetX());
        mcPrimaryVtxY.push_back(pMCPrimary->GetVertex().GetY());
        mcPrimaryVtxZ.push_back(pMCPrimary->GetVertex().GetZ());
        mcPrimaryEndX.push_back(pMCPrimary->GetEndpoint().GetX());
        mcPrimaryEndY.push_back(pMCPrimary->GetEndpoint().GetY());
        mcPrimaryEndZ.push_back(pMCPrimary->GetEndpoint().GetZ());
        nMCHitsTotal.push_back(mcPrimaryHitList.size());
        nMCHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList));
        nMCHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList));
        nMCHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList));

        int matchIndex(0), nPrimaryMatches(0), nPrimaryNuMatches(0), nPrimaryCRMatches(0), nPrimaryGoodNuMatches(0), nPrimaryNuSplits(0);
#ifdef MONITORING
        float recoVertexX(std::numeric_limits<float>::max()), recoVertexY(std::numeric_limits<float>::max()), recoVertexZ(std::numeric_limits<float>::max());
#endif
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));
            const bool isRecoTestBeam(LArPfoHelper::IsTestBeam(pfoToSharedHits.first));
            const bool isRecoTestBeamHierarchy(LArPfoHelper::IsTestBeam(LArPfoHelper::GetParentPfo(pfoToSharedHits.first)));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            // Tier (0) : Primary, (1) : Daughter, (2) : Granddaughter etc...  Note that the tier only increases for visible particle
            const int pfoHierarchyTier(LArPfoHelper::GetHierarchyTier(pfoToSharedHits.first));
            const int pfoId(pfoToIdMap.at(pfoToSharedHits.first));
            const int recoNuId(isRecoNeutrinoFinalState ? neutrinoPfoToIdMap.at(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first)) : -1);
            const int recoTBId(isRecoTestBeam || isRecoTestBeamHierarchy ? testBeamPfoToIdMap.at(LArPfoHelper::GetParentPfo(pfoToSharedHits.first)) : -1);

            if (0 == matchIndex++)
            {
                bestMatchPfoId.push_back(pfoId);
                bestMatchPfoPdg.push_back(pfoToSharedHits.first->GetParticleId());
                bestMatchPfoTier.push_back(pfoHierarchyTier);
                bestMatchPfoIsRecoNu.push_back(isRecoNeutrinoFinalState ? 1 : 0);
                bestMatchPfoRecoNuId.push_back(recoNuId);
                bestMatchPfoIsTestBeam.push_back(isRecoTestBeam ? 1 : 0);
                bestMatchPfoRecoTBId.push_back(recoTBId);
                bestMatchPfoNHitsTotal.push_back(pfoHitList.size());
                bestMatchPfoNHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                bestMatchPfoNHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                bestMatchPfoNHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                bestMatchPfoNSharedHitsTotal.push_back(sharedHitList.size());
                bestMatchPfoNSharedHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                bestMatchPfoNSharedHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                bestMatchPfoNSharedHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));
#ifdef MONITORING
                try
                {
                    const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(isRecoNeutrinoFinalState ? LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first) : pfoToSharedHits.first));
                    recoVertexX = pRecoVertex->GetPosition().GetX();
                    recoVertexY = pRecoVertex->GetPosition().GetY();
                    recoVertexZ = pRecoVertex->GetPosition().GetZ();
                }
                catch (const StatusCodeException &) {}
#endif
            }

            if (isGoodMatch) ++nPrimaryMatches;

            if (isRecoNeutrinoFinalState)
            {
                const Pfo *const pRecoNeutrino(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first));
                const bool isSplitRecoNeutrino(!recoNeutrinos.empty() && !recoNeutrinos.count(pRecoNeutrino));
                if (!isSplitRecoNeutrino && isGoodMatch) ++nPrimaryGoodNuMatches;
                if (isSplitRecoNeutrino && isBeamNeutrinoFinalState && isGoodMatch) ++nPrimaryNuSplits;
                recoNeutrinos.insert(pRecoNeutrino);
            }

            if (!m_testBeamMode && !m_testBeamHierarchyMode)
            {
                // If not in test beam mode proceed as standard
                if (isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryNuMatches;
                if (!isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryCRMatches;
            }
            else if (m_testBeamHierarchyMode)
            {
                // ATTN: In hierarchy mode let NuMatches become effective TBHierarchyMatches and treat the same
                if (isRecoTestBeamHierarchy && isGoodMatch) ++nPrimaryNuMatches;
                if (!isRecoTestBeamHierarchy && isGoodMatch) ++nPrimaryCRMatches;

                // Account for splitting of test beam particle into separate reconstructed primary pfos
                const Pfo *const pRecoTB(LArPfoHelper::GetParentPfo(pfoToSharedHits.first));
                const bool isSplitRecoTBHierarchy(!recoTestBeamHierarchies.empty() && !recoTestBeamHierarchies.count(pRecoTB));
                if (!isSplitRecoTBHierarchy && isGoodMatch) ++nPrimaryGoodNuMatches;
                if (isSplitRecoTBHierarchy && isLeadingBeamParticle && isGoodMatch) ++nPrimaryNuSplits;
                recoTestBeamHierarchies.insert(pRecoTB);
            }
            else
            {
                if (isRecoTestBeam && isGoodMatch) ++nPrimaryNuMatches;
                if (!isRecoTestBeam && isGoodMatch) ++nPrimaryCRMatches;
            }

            if (m_testBeamHierarchyMode) for (int tier = 0; tier < mcHierarchyTier; tier++) targetSS << "    ";

            targetSS << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                     << "MatchedPfoId " << pfoId
                     << ", Nu " << isRecoNeutrinoFinalState;
            if (isRecoNeutrinoFinalState) targetSS << " [NuId: " << recoNuId << "]";
            targetSS << ", TB " << isRecoTestBeam;
            if (m_testBeamHierarchyMode) targetSS << ", TB Hierarchy " << isRecoTestBeamHierarchy;
            if (isRecoTestBeamHierarchy) targetSS << " [TBId: " << recoTBId << "]";
            targetSS << ", CR " << (!isRecoNeutrinoFinalState && !isRecoTestBeam && !isRecoTestBeamHierarchy)
                     << ", PDG " << pfoToSharedHits.first->GetParticleId();
            if (m_testBeamHierarchyMode) targetSS << ", Tier " << pfoHierarchyTier;
            targetSS << ", nMatchedHits " << sharedHitList.size()
                     << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                     << ", nPfoHits " << pfoHitList.size()
                     << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;
        }

        if (mcToPfoHitSharingMap.at(pMCPrimary).empty())
        {
            if (m_testBeamHierarchyMode) for (int tier = 0; tier < mcHierarchyTier; tier++) targetSS << "    ";
            targetSS << "-No matched Pfo" << std::endl;
            bestMatchPfoId.push_back(-1); bestMatchPfoPdg.push_back(0); bestMatchPfoTier.push_back(-1); bestMatchPfoIsRecoNu.push_back(0); bestMatchPfoRecoNuId.push_back(-1); bestMatchPfoIsTestBeam.push_back(0); bestMatchPfoRecoTBId.push_back(-1);
            bestMatchPfoNHitsTotal.push_back(0); bestMatchPfoNHitsU.push_back(0); bestMatchPfoNHitsV.push_back(0); bestMatchPfoNHitsW.push_back(0);
            bestMatchPfoNSharedHitsTotal.push_back(0); bestMatchPfoNSharedHitsU.push_back(0); bestMatchPfoNSharedHitsV.push_back(0); bestMatchPfoNSharedHitsW.push_back(0);
        }

        nPrimaryMatchedPfos.push_back(nPrimaryMatches);
        nPrimaryMatchedNuPfos.push_back(nPrimaryNuMatches);
        nPrimaryMatchedCRPfos.push_back(nPrimaryCRMatches);
        nTargetMatches += nPrimaryMatches;
        nTargetNuMatches += nPrimaryNuMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodNuMatches += nPrimaryGoodNuMatches;
        nTargetNuSplits += nPrimaryNuSplits;
        if (0 == nPrimaryMatches) ++nTargetNuLosses;

	if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNuanceCode", mcNuanceCode));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isNeutrino", isBeamNeutrinoFinalState));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isBeamParticle", isBeamParticle));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCosmicRay", isCosmicRay));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetPrimaries", nTargetPrimaries));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexX", targetVertexX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexY", targetVertexY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexZ", targetVertexZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexX", recoVertexX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexY", recoVertexY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexZ", recoVertexZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryId", &mcPrimaryId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPdg", &mcPrimaryPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryTier", &mcPrimaryTier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryE", &mcPrimaryE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPX", &mcPrimaryPX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPY", &mcPrimaryPY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPZ", &mcPrimaryPZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxX", &mcPrimaryVtxX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxY", &mcPrimaryVtxY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxZ", &mcPrimaryVtxZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndX", &mcPrimaryEndX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndY", &mcPrimaryEndY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndZ", &mcPrimaryEndZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsTotal", &nMCHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsU", &nMCHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsV", &nMCHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsW", &nMCHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedPfos", &nPrimaryMatchedPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedNuPfos", &nPrimaryMatchedNuPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedCRPfos", &nPrimaryMatchedCRPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoId", &bestMatchPfoId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPdg", &bestMatchPfoPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTier", &bestMatchPfoTier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsRecoNu", &bestMatchPfoIsRecoNu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoRecoNuId", &bestMatchPfoRecoNuId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsTotal", &bestMatchPfoNHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsU", &bestMatchPfoNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsV", &bestMatchPfoNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsW", &bestMatchPfoNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsTotal", &bestMatchPfoNSharedHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsU", &bestMatchPfoNSharedHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsV", &bestMatchPfoNSharedHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsW", &bestMatchPfoNSharedHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetMatches", nTargetMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuMatches", nTargetNuMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetCRMatches", nTargetCRMatches));

            if (m_testBeamMode || m_testBeamHierarchyMode)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsTestBeam", &bestMatchPfoIsTestBeam));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoRecoTBId", &bestMatchPfoRecoTBId));
            }

            if (!m_testBeamMode || m_testBeamHierarchyMode)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetGoodNuMatches", nTargetGoodNuMatches));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuSplits", nTargetNuSplits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuLosses", nTargetNuLosses));
            }
        }

        if ((isLastNeutrinoPrimary && !m_testBeamHierarchyMode) || (isBeamParticle && m_testBeamMode) || isCosmicRay || (isLastTestBeamLeading && m_testBeamHierarchyMode))
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries, m_testBeamHierarchyMode));
#ifdef MONITORING
            const int interactionTypeInt(static_cast<int>(interactionType));
#endif
            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectNu(isBeamNeutrinoFinalState && (nTargetGoodNuMatches == nTargetNuMatches) && (nTargetGoodNuMatches == nTargetPrimaries) && (nTargetCRMatches == 0) && (nTargetNuSplits == 0) && (nTargetNuLosses == 0));
            const int isCorrectTB(isBeamParticle && (nTargetNuMatches == 1) && (nTargetCRMatches == 0));
            const int isCorrectTBHierarchy(isLeadingBeamParticle && (nTargetGoodNuMatches == nTargetNuMatches) && (nTargetGoodNuMatches == nTargetPrimaries) && (nTargetCRMatches == 0) && (nTargetNuSplits == 0) && (nTargetNuLosses == 0));
            const int isCorrectCR(isCosmicRay && (nTargetNuMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeNu(isCosmicRay && (nTargetNuMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitNu(!isCosmicRay && ((nTargetNuMatches > nTargetPrimaries) || (nTargetNuSplits > 0)));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            std::stringstream outcomeSS;
            outcomeSS << LArInteractionTypeHelper::ToString(interactionType) << " (Nuance " << mcNuanceCode << ", Nu " << isBeamNeutrinoFinalState << ", TB " << isBeamParticle << ", CR " << isCosmicRay << ")" << std::endl;

            if (isLastNeutrinoPrimary) ++nTotalNu;
            if (isBeamParticle) ++nTotalTB;
            if (isLastTestBeamLeading) ++nTotalTBHierarchy;
            if (isCosmicRay) ++nTotalCR;
            if (isCorrectNu) ++nCorrectNu;
            if (isCorrectTB) ++nCorrectTB;
            if (isCorrectTBHierarchy) ++nCorrectTBHierarchy;
            if (isCorrectCR) ++nCorrectCR;
            if (isFakeNu) ++nFakeNu;
            if (isFakeCR) ++nFakeCR;
            if (isSplitNu) ++nSplitNu;
            if (isSplitCR) ++nSplitCR;
            if (isLost) ++nLost;

            if (isCorrectNu && !m_testBeamHierarchyMode) outcomeSS << "IsCorrectNu ";
            if (isCorrectTB && !m_testBeamHierarchyMode) outcomeSS << "IsCorrectTB ";
            if (isCorrectTBHierarchy && m_testBeamHierarchyMode) outcomeSS << "IsCorrectTBHierarchy";
            if (isCorrectCR) outcomeSS << "IsCorrectCR ";
            if (isFakeNu) outcomeSS << "IsFake" << name << " ";
            if (isFakeCR) outcomeSS << "IsFakeCR ";
            if (isSplitNu) outcomeSS << "isSplit" << name << " ";
            if (isSplitCR) outcomeSS << "IsSplitCR ";
            if (isLost) outcomeSS << "IsLost ";
            if (nTargetNuMatches > 0) outcomeSS << "(N" << name << "Matches: " << nTargetNuMatches << ") ";
            if (nTargetNuLosses > 0) outcomeSS << "(N" << name << "Losses: " << nTargetNuLosses << ") ";
            if (nTargetNuSplits > 0) outcomeSS << "(N" << name << "Splits: " << nTargetNuSplits << ") ";
            if (nTargetCRMatches > 0) outcomeSS << "(NCRMatches: " << nTargetCRMatches << ") ";
            if (printToScreen) std::cout << outcomeSS.str() << std::endl << targetSS.str() << std::endl;

            if (fillTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionTypeInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectNu", isCorrectNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectTB", isCorrectTB));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectTBHierarchy", isCorrectTBHierarchy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectCR", isCorrectCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeNu", isFakeNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeCR", isFakeCR));

                if (!m_testBeamMode || m_testBeamHierarchyMode)
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitNu", isSplitNu));

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitCR", isSplitCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isLost", isLost));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }

            targetSS.str(std::string()); targetSS.clear();
            recoNeutrinos.clear(); associatedMCPrimaries.clear();
            nTargetMatches = 0; nTargetNuMatches = 0; nTargetCRMatches = 0; nTargetGoodNuMatches = 0; nTargetNuSplits = 0; nTargetNuLosses = 0;
            mcPrimaryId.clear(); mcPrimaryPdg.clear(); mcPrimaryTier.clear(); nMCHitsTotal.clear(); nMCHitsU.clear(); nMCHitsV.clear(); nMCHitsW.clear();
            mcPrimaryE.clear(); mcPrimaryPX.clear(); mcPrimaryPY.clear(); mcPrimaryPZ.clear();
            mcPrimaryVtxX.clear(); mcPrimaryVtxY.clear(); mcPrimaryVtxZ.clear(); mcPrimaryEndX.clear(); mcPrimaryEndY.clear(); mcPrimaryEndZ.clear();
            nPrimaryMatchedPfos.clear(); nPrimaryMatchedNuPfos.clear(); nPrimaryMatchedCRPfos.clear();
            bestMatchPfoId.clear(); bestMatchPfoPdg.clear(); bestMatchPfoTier.clear(); bestMatchPfoIsRecoNu.clear(); bestMatchPfoRecoNuId.clear(); bestMatchPfoIsTestBeam.clear(); bestMatchPfoRecoTBId.clear();
            bestMatchPfoNHitsTotal.clear(); bestMatchPfoNHitsU.clear(); bestMatchPfoNHitsV.clear(); bestMatchPfoNHitsW.clear();
            bestMatchPfoNSharedHitsTotal.clear(); bestMatchPfoNSharedHitsU.clear(); bestMatchPfoNSharedHitsV.clear(); bestMatchPfoNSharedHitsW.clear();
        }
    }

    if (useInterpretedMatching)
    {
        std::stringstream summarySS;
        summarySS << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalNu > 0 && !m_testBeamHierarchyMode) summarySS << "#CorrectNu: " << nCorrectNu << "/" << nTotalNu << ", Fraction: " << (nTotalNu > 0 ? static_cast<float>(nCorrectNu) / static_cast<float>(nTotalNu) : 0.f) << std::endl;
        if (nTotalTB > 0  && !m_testBeamHierarchyMode) summarySS << "#CorrectTB: " << nCorrectTB << "/" << nTotalTB << ", Fraction: " << (nTotalTB > 0 ? static_cast<float>(nCorrectTB) / static_cast<float>(nTotalTB) : 0.f) << std::endl;
        if (nTotalTBHierarchy > 0 && m_testBeamHierarchyMode) summarySS << "#CorrectTBHierarchy: " << nCorrectTBHierarchy << "/" << nTotalTBHierarchy << ", Fraction: " << (nTotalTBHierarchy > 0 ? static_cast<float>(nCorrectTBHierarchy) / static_cast<float>(nTotalTBHierarchy) : 0.f) << std::endl;
        if (nTotalCR > 0) summarySS << "#CorrectCR: " << nCorrectCR << "/" << nTotalCR << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if (nFakeNu > 0) summarySS << "#Fake" << name << ": " << nFakeNu << " ";
        if (nFakeCR > 0) summarySS << "#FakeCR: " << nFakeCR << " ";
        if (nSplitNu > 0) summarySS << "#Split" << name << ": " << nSplitNu << " ";
        if (nSplitCR > 0) summarySS << "#SplitCR: " << nSplitCR << " ";
        if (nLost > 0) summarySS << "#Lost: " << nLost << " ";
        if (nFakeNu || nFakeCR || nSplitNu || nSplitCR || nLost) summarySS << std::endl;
        if (printToScreen) std::cout << summarySS.str();
    }

    if (printToScreen) std::cout << "------------------------------------------------------------------------------------------------" << std::endl << std::endl;
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

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
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

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
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
        "TestBeamMode", m_testBeamMode));

     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TestBeamHierarchyMode", m_testBeamHierarchyMode));

    if (m_testBeamMode && m_testBeamHierarchyMode)
    {
        std::cout << "EventValidationAlgorithm::ReadSettings - TestBeamMode or TestBeamHierarchyMode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

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
