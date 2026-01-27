/**
 *  @file   larpandoracontent/LArMonitoring/NeutrinoEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/NeutrinoEventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

NeutrinoEventValidationAlgorithm::NeutrinoEventValidationAlgorithm() :
    m_useTrueNeutrinosOnly(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoEventValidationAlgorithm::~NeutrinoEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList,
    const CaloHitList *const pCaloHitList, const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly)
            LArMCParticleHelper::SelectReconstructableMCParticles(
                pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        LArMCParticleHelper::PrimaryParameters parameters(m_primaryParameters);
        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly)
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
            if (LArPfoHelper::IsFinalState(pPfo))
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

    // ATTN : Ensure all mc primaries have an entry in mcToPfoHitSharingMap, even if no associated pfos.
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetAllMCParticleToHitsMap()}, mcPrimaryVector);
    for (const MCParticle *pMCParticle : mcPrimaryVector)
    {
        if (mcToPfoHitSharingMap.find(pMCParticle) == mcToPfoHitSharingMap.end())
        {
            LArMCParticleHelper::PfoToSharedHitsVector pfoToSharedHitsVector;
            mcToPfoHitSharingMap.insert(LArMCParticleHelper::MCParticleToPfoHitSharingMap::value_type(pMCParticle, pfoToSharedHitsVector));
        }
    }

    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoEventValidationAlgorithm::ProcessOutput(
    const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    if (printToScreen && useInterpretedMatching)
        std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen)
        std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(
        useInterpretedMatching ? validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    // Neutrino Validation Bookkeeping
    int nNeutrinoPrimaries(0);
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            ++nNeutrinoPrimaries;

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    int pfoIndex(0), neutrinoPfoIndex(0);
    PfoToIdMap pfoToIdMap, neutrinoPfoToIdMap;

    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));
        const Pfo *const pRecoNeutrino(LArPfoHelper::IsNeutrinoFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentNeutrino(pPrimaryPfo) : nullptr);

        if (pRecoNeutrino && !neutrinoPfoToIdMap.count(pRecoNeutrino))
            neutrinoPfoToIdMap.insert(PfoToIdMap::value_type(pRecoNeutrino, ++neutrinoPfoIndex));
    }

    LArMCParticleHelper::MCParticleIntMap triggeredToLeading, triggeredToLeadingCounter;

    PfoSet recoNeutrinos;
    MCParticleList associatedMCPrimaries;

    int nCorrectNu(0), nTotalNu(0), nCorrectCR(0), nTotalCR(0);
    int nFakeNu(0), nFakeCR(0), nSplitNu(0), nSplitCR(0), nLost(0), mcPrimaryIndex(0), nTargetMatches(0), nTargetNuMatches(0);
    int nTargetCRMatches(0), nTargetGoodNuMatches(0), nTargetNuSplits(0), nTargetNuLosses(0);
    IntVector mcPrimaryId, mcPrimaryPdg, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedNuPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoIsRecoNu, bestMatchPfoRecoNuId;
    IntVector bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;

    std::stringstream targetSS;
    const std::string name("Nu");

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!isTargetPrimary && !hasMatch)
            continue;

        associatedMCPrimaries.push_back(pMCPrimary);
        const int nTargetPrimaries(associatedMCPrimaries.size());
        const bool isLastNeutrinoPrimary(++mcPrimaryIndex == nNeutrinoPrimaries);
        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        const int isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));
#ifdef MONITORING
        const CartesianVector &targetVertex(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetVertex());
        const float targetVertexX(targetVertex.GetX()), targetVertexY(targetVertex.GetY()), targetVertexZ(targetVertex.GetZ());
#endif

        targetSS << (!isTargetPrimary ? "(Non target) " : "") << "PrimaryId " << mcPrimaryIndex << ", Nu " << isBeamNeutrinoFinalState
                 << ", CR " << isCosmicRay << ", MCPDG " << pMCPrimary->GetParticleId() << ", Energy " << pMCPrimary->GetEnergy()
                 << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude() << ", nMCHits "
                 << mcPrimaryHitList.size() << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList) << ", "
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

        int matchIndex(0), nPrimaryMatches(0), nPrimaryNuMatches(0), nPrimaryCRMatches(0), nPrimaryGoodNuMatches(0), nPrimaryNuSplits(0);
#ifdef MONITORING
        float recoVertexX(std::numeric_limits<float>::max()), recoVertexY(std::numeric_limits<float>::max()),
            recoVertexZ(std::numeric_limits<float>::max());
#endif
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            const int pfoId(pfoToIdMap.at(pfoToSharedHits.first));
            const int recoNuId(
                isRecoNeutrinoFinalState ? static_cast<int>(neutrinoPfoToIdMap.at(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first))) : -1);

            if (0 == matchIndex++)
            {
                bestMatchPfoId.push_back(pfoId);
                bestMatchPfoPdg.push_back(pfoToSharedHits.first->GetParticleId());
                bestMatchPfoIsRecoNu.push_back(isRecoNeutrinoFinalState ? 1 : 0);
                bestMatchPfoRecoNuId.push_back(recoNuId);
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
                    const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(
                        isRecoNeutrinoFinalState ? LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first) : pfoToSharedHits.first));
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

            if (isRecoNeutrinoFinalState)
            {
                const Pfo *const pRecoNeutrino(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first));
                const bool isSplitRecoNeutrino(!recoNeutrinos.empty() && !recoNeutrinos.count(pRecoNeutrino));
                if (!isSplitRecoNeutrino && isGoodMatch)
                    ++nPrimaryGoodNuMatches;
                if (isSplitRecoNeutrino && isBeamNeutrinoFinalState && isGoodMatch)
                    ++nPrimaryNuSplits;
                recoNeutrinos.insert(pRecoNeutrino);
            }

            if (isRecoNeutrinoFinalState && isGoodMatch)
                ++nPrimaryNuMatches;
            if (!isRecoNeutrinoFinalState && isGoodMatch)
                ++nPrimaryCRMatches;

            targetSS << "-" << (!isGoodMatch ? "(Below threshold) " : "") << "MatchedPfoId " << pfoId << ", Nu " << isRecoNeutrinoFinalState;
            if (isRecoNeutrinoFinalState)
                targetSS << " [NuId: " << recoNuId << "]";
            targetSS << ", CR " << (!isRecoNeutrinoFinalState) << ", PDG " << pfoToSharedHits.first->GetParticleId() << ", nMatchedHits "
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
            bestMatchPfoIsRecoNu.push_back(0);
            bestMatchPfoRecoNuId.push_back(-1);
            bestMatchPfoNHitsTotal.push_back(0);
            bestMatchPfoNHitsU.push_back(0);
            bestMatchPfoNHitsV.push_back(0);
            bestMatchPfoNHitsW.push_back(0);
            bestMatchPfoNSharedHitsTotal.push_back(0);
            bestMatchPfoNSharedHitsU.push_back(0);
            bestMatchPfoNSharedHitsV.push_back(0);
            bestMatchPfoNSharedHitsW.push_back(0);
        }

        nPrimaryMatchedPfos.push_back(nPrimaryMatches);
        nPrimaryMatchedNuPfos.push_back(nPrimaryNuMatches);
        nPrimaryMatchedCRPfos.push_back(nPrimaryCRMatches);
        nTargetMatches += nPrimaryMatches;
        nTargetNuMatches += nPrimaryNuMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodNuMatches += nPrimaryGoodNuMatches;
        nTargetNuSplits += nPrimaryNuSplits;
        if (0 == nPrimaryMatches)
            ++nTargetNuLosses;

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNuanceCode", mcNuanceCode));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isNeutrino", isBeamNeutrinoFinalState));
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedNuPfos", &nPrimaryMatchedNuPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedCRPfos", &nPrimaryMatchedCRPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoId", &bestMatchPfoId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPdg", &bestMatchPfoPdg));
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetGoodNuMatches", nTargetGoodNuMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuSplits", nTargetNuSplits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuLosses", nTargetNuLosses));
        }

        if (isLastNeutrinoPrimary || isCosmicRay)
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries));
#ifdef MONITORING
            const int interactionTypeInt(static_cast<int>(interactionType));
#endif
            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectNu(isBeamNeutrinoFinalState && (nTargetGoodNuMatches == nTargetNuMatches) &&
                (nTargetGoodNuMatches == nTargetPrimaries) && (nTargetCRMatches == 0) && (nTargetNuSplits == 0) && (nTargetNuLosses == 0));
            const int isCorrectCR(isCosmicRay && (nTargetNuMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeNu(isCosmicRay && (nTargetNuMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitNu(!isCosmicRay && ((nTargetNuMatches > nTargetPrimaries) || (nTargetNuSplits > 0)));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            std::stringstream outcomeSS;
            outcomeSS << LArInteractionTypeHelper::ToString(interactionType) << " (Nuance " << mcNuanceCode << ", Nu "
                      << isBeamNeutrinoFinalState << ", CR " << isCosmicRay << ")" << std::endl;

            if (isLastNeutrinoPrimary)
                ++nTotalNu;
            if (isCosmicRay)
                ++nTotalCR;
            if (isCorrectNu)
                ++nCorrectNu;
            if (isCorrectCR)
                ++nCorrectCR;
            if (isFakeNu)
                ++nFakeNu;
            if (isFakeCR)
                ++nFakeCR;
            if (isSplitNu)
                ++nSplitNu;
            if (isSplitCR)
                ++nSplitCR;
            if (isLost)
                ++nLost;

            if (isCorrectNu)
                outcomeSS << "IsCorrectNu ";
            if (isCorrectCR)
                outcomeSS << "IsCorrectCR ";
            if (isFakeNu)
                outcomeSS << "IsFake" << name << " ";
            if (isFakeCR)
                outcomeSS << "IsFakeCR ";
            if (isSplitNu)
                outcomeSS << "isSplit" << name << " ";
            if (isSplitCR)
                outcomeSS << "IsSplitCR ";
            if (isLost)
                outcomeSS << "IsLost ";
            if (nTargetNuMatches > 0)
                outcomeSS << "(N" << name << "Matches: " << nTargetNuMatches << ") ";
            if (nTargetNuLosses > 0)
                outcomeSS << "(N" << name << "Losses: " << nTargetNuLosses << ") ";
            if (nTargetNuSplits > 0)
                outcomeSS << "(N" << name << "Splits: " << nTargetNuSplits << ") ";
            if (nTargetCRMatches > 0)
                outcomeSS << "(NCRMatches: " << nTargetCRMatches << ") ";
            if (printToScreen)
                std::cout << outcomeSS.str() << std::endl << targetSS.str() << std::endl;

            if (fillTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionTypeInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectNu", isCorrectNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectCR", isCorrectCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeNu", isFakeNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeCR", isFakeCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitNu", isSplitNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitCR", isSplitCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isLost", isLost));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }

            targetSS.str(std::string());
            targetSS.clear();
            recoNeutrinos.clear();
            associatedMCPrimaries.clear();
            nTargetMatches = 0;
            nTargetNuMatches = 0;
            nTargetCRMatches = 0;
            nTargetGoodNuMatches = 0;
            nTargetNuSplits = 0;
            nTargetNuLosses = 0;
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
            nPrimaryMatchedNuPfos.clear();
            nPrimaryMatchedCRPfos.clear();
            bestMatchPfoId.clear();
            bestMatchPfoPdg.clear();
            bestMatchPfoIsRecoNu.clear();
            bestMatchPfoRecoNuId.clear();
            bestMatchPfoNHitsTotal.clear();
            bestMatchPfoNHitsU.clear();
            bestMatchPfoNHitsV.clear();
            bestMatchPfoNHitsW.clear();
            bestMatchPfoNSharedHitsTotal.clear();
            bestMatchPfoNSharedHitsU.clear();
            bestMatchPfoNSharedHitsV.clear();
            bestMatchPfoNSharedHitsW.clear();
        }
    }

    if (useInterpretedMatching)
    {
        std::stringstream summarySS;
        summarySS << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalNu > 0)
            summarySS << "#CorrectNu: " << nCorrectNu << "/" << nTotalNu
                      << ", Fraction: " << (nTotalNu > 0 ? static_cast<float>(nCorrectNu) / static_cast<float>(nTotalNu) : 0.f) << std::endl;
        if (nTotalCR > 0)
            summarySS << "#CorrectCR: " << nCorrectCR << "/" << nTotalCR
                      << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if (nFakeNu > 0)
            summarySS << "#Fake" << name << ": " << nFakeNu << " ";
        if (nFakeCR > 0)
            summarySS << "#FakeCR: " << nFakeCR << " ";
        if (nSplitNu > 0)
            summarySS << "#Split" << name << ": " << nSplitNu << " ";
        if (nSplitCR > 0)
            summarySS << "#SplitCR: " << nSplitCR << " ";
        if (nLost > 0)
            summarySS << "#Lost: " << nLost << " ";
        if (nFakeNu || nFakeCR || nSplitNu || nSplitCR || nLost)
            summarySS << std::endl;
        if (printToScreen)
            std::cout << summarySS.str();
    }

    if (printToScreen)
        std::cout << "------------------------------------------------------------------------------------------------" << std::endl
                  << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
