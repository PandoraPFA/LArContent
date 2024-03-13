/**
 *  @file   larpandoracontent/LArMonitoring/TestBeamEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the test beam event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/TestBeamEventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

TestBeamEventValidationAlgorithm::TestBeamEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TestBeamEventValidationAlgorithm::~TestBeamEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList,
    const CaloHitList *const pCaloHitList, const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        LArMCParticleHelper::PrimaryParameters parameters(m_primaryParameters);
        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, allMCParticleToHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

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
            if (pPfo->GetParentPfoList().empty())
                finalStatePfos.push_back(pPfo);
        }

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(
            finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, m_primaryParameters.m_foldBackHierarchy);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(
        validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamEventValidationAlgorithm::ProcessOutput(
    const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    static int eventNumber{-1};
    ++eventNumber;
    if (printToScreen && useInterpretedMatching)
        std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen)
        std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(
        useInterpretedMatching ? validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    // Test Beam Validation Bookkeeping
    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    int pfoIndex(0), testBeamPfoIndex(0);
    PfoToIdMap pfoToIdMap, testBeamPfoToIdMap;

    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));
        const Pfo *const pRecoTestBeam(LArPfoHelper::IsTestBeamFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentPfo(pPrimaryPfo) : nullptr);

        if (pRecoTestBeam && !testBeamPfoToIdMap.count(pRecoTestBeam))
            testBeamPfoToIdMap.insert(PfoToIdMap::value_type(pRecoTestBeam, ++testBeamPfoIndex));
    }

    LArMCParticleHelper::MCParticleIntMap triggeredToLeading, triggeredToLeadingCounter;

    MCParticleList associatedMCPrimaries;

    int nCorrectTB(0), nTotalTB(0), nCorrectCR(0), nTotalCR(0), nFakeTB(0), nFakeCR(0), nSplitTB(0), nSplitCR(0), nLost(0);
    int mcPrimaryIndex(0), nTargetMatches(0), nTargetTBMatches(0), nTargetCRMatches(0), nTargetGoodTBMatches(0);
    IntVector mcPrimaryId, mcPrimaryPdg, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedTBPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoIsTB;
    IntVector bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;
    FloatVector bestMatchPfoX0;

    std::stringstream targetSS;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!isTargetPrimary && !hasMatch)
            continue;

        associatedMCPrimaries.push_back(pMCPrimary);
        ++mcPrimaryIndex;
        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        const int isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));
#ifdef MONITORING
        const int nTargetPrimaries(associatedMCPrimaries.size());
        const CartesianVector &targetVertex(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetVertex());
        const float targetVertexX(targetVertex.GetX()), targetVertexY(targetVertex.GetY()), targetVertexZ(targetVertex.GetZ());
#endif

        targetSS << (!isTargetPrimary ? "(Non target) " : "") << "PrimaryId " << mcPrimaryIndex << ", TB " << isBeamParticle << ", CR "
                 << isCosmicRay << ", MCPDG " << pMCPrimary->GetParticleId() << ", Energy " << pMCPrimary->GetEnergy() << ", Dist. "
                 << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude() << ", nMCHits " << mcPrimaryHitList.size() << " ("
                 << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList) << ", "
                 << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList) << ", "
                 << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        mcPrimaryId.push_back(mcPrimaryIndex);
        mcPrimaryPdg.push_back(pMCPrimary->GetParticleId());
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

        int matchIndex(0), nPrimaryMatches(0), nPrimaryTBMatches(0), nPrimaryCRMatches(0), nPrimaryGoodNuMatches(0);
#ifdef MONITORING
        float recoVertexX(std::numeric_limits<float>::max()), recoVertexY(std::numeric_limits<float>::max()),
            recoVertexZ(std::numeric_limits<float>::max());
#endif
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoTestBeam(LArPfoHelper::IsTestBeam(pfoToSharedHits.first));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            const int pfoId(pfoToIdMap.at(pfoToSharedHits.first));

            if (0 == matchIndex++)
            {
                bestMatchPfoId.push_back(pfoId);
                bestMatchPfoPdg.push_back(pfoToSharedHits.first->GetParticleId());
                bestMatchPfoIsTB.push_back(isRecoTestBeam ? 1 : 0);
                bestMatchPfoNHitsTotal.push_back(pfoHitList.size());
                bestMatchPfoNHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                bestMatchPfoNHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                bestMatchPfoNHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                bestMatchPfoNSharedHitsTotal.push_back(sharedHitList.size());
                bestMatchPfoNSharedHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                bestMatchPfoNSharedHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                bestMatchPfoNSharedHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));
                bestMatchPfoX0.push_back(pfoToSharedHits.first->GetPropertiesMap().count("X0") ? pfoToSharedHits.first->GetPropertiesMap().at("X0")
                                                                                               : std::numeric_limits<float>::max());
#ifdef MONITORING
                try
                {
                    const Vertex *const pRecoVertex(isRecoTestBeam ? LArPfoHelper::GetTestBeamInteractionVertex(pfoToSharedHits.first)
                                                                   : LArPfoHelper::GetVertex(pfoToSharedHits.first));
                    recoVertexX = pRecoVertex->GetPosition().GetX();
                    recoVertexY = pRecoVertex->GetPosition().GetY();
                    recoVertexZ = pRecoVertex->GetPosition().GetZ();
                }
                catch (const StatusCodeException &)
                {
                }
#endif
            }

            if (isGoodMatch)
                ++nPrimaryMatches;

            if (isRecoTestBeam && isGoodMatch)
                ++nPrimaryTBMatches;
            if (!isRecoTestBeam && isGoodMatch)
                ++nPrimaryCRMatches;

            targetSS << "-" << (!isGoodMatch ? "(Below threshold) " : "") << "MatchedPfoId " << pfoId << ", TB " << isRecoTestBeam
                     << ", CR " << (!isRecoTestBeam) << ", PDG " << pfoToSharedHits.first->GetParticleId() << ", nMatchedHits "
                     << sharedHitList.size() << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList) << ", "
                     << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList) << ", "
                     << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                     << ", nPfoHits " << pfoHitList.size() << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList) << ", "
                     << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList) << ", "
                     << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;
        }

        if (mcToPfoHitSharingMap.at(pMCPrimary).empty())
        {
            targetSS << "-No matched Pfo" << std::endl;
            bestMatchPfoId.push_back(-1);
            bestMatchPfoPdg.push_back(0);
            bestMatchPfoIsTB.push_back(0);
            bestMatchPfoNHitsTotal.push_back(0);
            bestMatchPfoNHitsU.push_back(0);
            bestMatchPfoNHitsV.push_back(0);
            bestMatchPfoNHitsW.push_back(0);
            bestMatchPfoNSharedHitsTotal.push_back(0);
            bestMatchPfoNSharedHitsU.push_back(0);
            bestMatchPfoNSharedHitsV.push_back(0);
            bestMatchPfoNSharedHitsW.push_back(0);
            bestMatchPfoX0.push_back(std::numeric_limits<float>::max());
        }

        nPrimaryMatchedPfos.push_back(nPrimaryMatches);
        nPrimaryMatchedTBPfos.push_back(nPrimaryTBMatches);
        nPrimaryMatchedCRPfos.push_back(nPrimaryCRMatches);
        nTargetMatches += nPrimaryMatches;
        nTargetTBMatches += nPrimaryTBMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodTBMatches += nPrimaryGoodNuMatches;

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", eventNumber));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNuanceCode", mcNuanceCode));
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedTBPfos", &nPrimaryMatchedTBPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedCRPfos", &nPrimaryMatchedCRPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoId", &bestMatchPfoId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPdg", &bestMatchPfoPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsTotal", &bestMatchPfoNHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsU", &bestMatchPfoNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsV", &bestMatchPfoNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsW", &bestMatchPfoNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsTotal", &bestMatchPfoNSharedHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsU", &bestMatchPfoNSharedHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsV", &bestMatchPfoNSharedHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsW", &bestMatchPfoNSharedHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoX0", &bestMatchPfoX0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetMatches", nTargetMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetTBMatches", nTargetTBMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetCRMatches", nTargetCRMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsTB", &bestMatchPfoIsTB));
        }

        if (isBeamParticle || isCosmicRay)
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries));
#ifdef MONITORING
            const int interactionTypeInt(static_cast<int>(interactionType));
#endif
            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectTB(isBeamParticle && (nTargetTBMatches == 1) && (nTargetCRMatches == 0));
            const int isCorrectCR(isCosmicRay && (nTargetTBMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeTB(isCosmicRay && (nTargetTBMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitTB(!isCosmicRay && (nTargetTBMatches > 1));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            std::stringstream outcomeSS;
            outcomeSS << LArInteractionTypeHelper::ToString(interactionType) << " (Nuance " << mcNuanceCode << ", TB " << isBeamParticle
                      << ", CR " << isCosmicRay << ")" << std::endl;

            if (isBeamParticle)
                ++nTotalTB;
            if (isCosmicRay)
                ++nTotalCR;
            if (isCorrectTB)
                ++nCorrectTB;
            if (isCorrectCR)
                ++nCorrectCR;
            if (isFakeTB)
                ++nFakeTB;
            if (isFakeCR)
                ++nFakeCR;
            if (isSplitTB)
                ++nSplitTB;
            if (isSplitCR)
                ++nSplitCR;
            if (isLost)
                ++nLost;

            if (isCorrectTB)
                outcomeSS << "IsCorrectTB ";
            if (isCorrectCR)
                outcomeSS << "IsCorrectCR ";
            if (isFakeTB)
                outcomeSS << "IsFakeTB ";
            if (isFakeCR)
                outcomeSS << "IsFakeCR ";
            if (isSplitTB)
                outcomeSS << "isSplitTB ";
            if (isSplitCR)
                outcomeSS << "IsSplitCR ";
            if (isLost)
                outcomeSS << "IsLost ";
            if (nTargetTBMatches > 0)
                outcomeSS << "(NTBMatches: " << nTargetTBMatches << ") ";
            if (nTargetCRMatches > 0)
                outcomeSS << "(NCRMatches: " << nTargetCRMatches << ") ";
            if (printToScreen)
                std::cout << outcomeSS.str() << std::endl << targetSS.str() << std::endl;

            if (fillTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionTypeInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectTB", isCorrectTB));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectCR", isCorrectCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeTB", isFakeTB));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeCR", isFakeCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitTB", isSplitTB));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitCR", isSplitCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isLost", isLost));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }

            targetSS.str(std::string());
            targetSS.clear();
            associatedMCPrimaries.clear();
            nTargetMatches = 0;
            nTargetTBMatches = 0;
            nTargetCRMatches = 0;
            nTargetGoodTBMatches = 0;
            mcPrimaryId.clear();
            mcPrimaryPdg.clear();
            nMCHitsTotal.clear();
            nMCHitsU.clear();
            nMCHitsV.clear();
            nMCHitsW.clear();
            mcPrimaryE.clear();
            mcPrimaryPX.clear();
            mcPrimaryPY.clear();
            mcPrimaryPZ.clear();
            mcPrimaryVtxX.clear();
            mcPrimaryVtxY.clear();
            mcPrimaryVtxZ.clear();
            mcPrimaryEndX.clear();
            mcPrimaryEndY.clear();
            mcPrimaryEndZ.clear();
            nPrimaryMatchedPfos.clear();
            nPrimaryMatchedTBPfos.clear();
            nPrimaryMatchedCRPfos.clear();
            bestMatchPfoId.clear();
            bestMatchPfoPdg.clear();
            bestMatchPfoIsTB.clear();
            bestMatchPfoNHitsTotal.clear();
            bestMatchPfoNHitsU.clear();
            bestMatchPfoNHitsV.clear();
            bestMatchPfoNHitsW.clear();
            bestMatchPfoNSharedHitsTotal.clear();
            bestMatchPfoNSharedHitsU.clear();
            bestMatchPfoNSharedHitsV.clear();
            bestMatchPfoNSharedHitsW.clear();
            bestMatchPfoX0.clear();
        }
    }

    if (useInterpretedMatching)
    {
        std::stringstream summarySS;
        summarySS << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalTB > 0)
            summarySS << "#CorrectTB: " << nCorrectTB << "/" << nTotalTB
                      << ", Fraction: " << (nTotalTB > 0 ? static_cast<float>(nCorrectTB) / static_cast<float>(nTotalTB) : 0.f) << std::endl;
        if (nTotalCR > 0)
            summarySS << "#CorrectCR: " << nCorrectCR << "/" << nTotalCR
                      << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if (nFakeTB > 0)
            summarySS << "#FakeTB: " << nFakeTB << " ";
        if (nFakeCR > 0)
            summarySS << "#FakeCR: " << nFakeCR << " ";
        if (nSplitTB > 0)
            summarySS << "#SplitTB: " << nSplitTB << " ";
        if (nSplitCR > 0)
            summarySS << "#SplitCR: " << nSplitCR << " ";
        if (nLost > 0)
            summarySS << "#Lost: " << nLost << " ";
        if (nFakeTB || nFakeCR || nSplitTB || nSplitCR || nLost)
            summarySS << std::endl;
        if (printToScreen)
            std::cout << summarySS.str();
    }

    if (printToScreen)
        std::cout << "------------------------------------------------------------------------------------------------" << std::endl
                  << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
