/**
 *  @file   larpandoracontent/LArMonitoring/TestBeamHierarchyEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the test beam hierarchy event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/TestBeamHierarchyEventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

TestBeamHierarchyEventValidationAlgorithm::TestBeamHierarchyEventValidationAlgorithm() :
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TestBeamHierarchyEventValidationAlgorithm::~TestBeamHierarchyEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamHierarchyEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::PrimaryParameters parameters;

        parameters.m_selectInputHits = m_selectInputHits;
        parameters.m_minHitSharingFraction = m_minHitSharingFraction;
        parameters.m_maxPhotonPropagation = m_maxPhotonPropagation;
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableTestBeamHierarchyMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsLeadingBeamParticle, targetMCParticleToHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableTestBeamHierarchyMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsLeadingBeamParticle, allMCParticleToHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

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

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap);
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

void TestBeamHierarchyEventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    // Test Beam Hierarchy Validation Pfo Bookkeeping
    int pfoIndex(0), testBeamPfoIndex(0);
    PfoToIdMap pfoToIdMap, testBeamPfoToIdMap;

    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));
        const Pfo *const pRecoTestBeam(LArPfoHelper::IsTestBeamFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentPfo(pPrimaryPfo) : nullptr); 

        if (pRecoTestBeam && !testBeamPfoToIdMap.count(pRecoTestBeam))
            testBeamPfoToIdMap.insert(PfoToIdMap::value_type(pRecoTestBeam, ++testBeamPfoIndex));
    }

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    // Test Beam Hierarchy Validation MCParticle Bookkeeping
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);
    LArMCParticleHelper::MCParticleIntMap triggeredToLeading, triggeredToLeadingCounter;

    // ATTN: At this stage the mcPrimaryVector is ordered from neutrinos, beam and then cosmics.  Here we extract the beam and reorder
    // to ensure the order follows primary parent beam 1, daughter 1 of beam 1, daughter 2 of beam 1, ..., primary parent beam 2,
    // daughter 1 of beam 2, etc... as expected by downstream logic
    MCParticleVector mcPrimaryVectorCopy(mcPrimaryVector), triggeredBeamParticles;
    LArMCParticleHelper::MCRelationMap leadingToTriggeredMap;
    mcPrimaryVector.clear();

    for (const MCParticle *const pMCPrimary : mcPrimaryVectorCopy)
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

    PfoSet recoTestBeamHierarchies;
    MCParticleList associatedMCPrimaries;

    int nCorrectTB(0), nTotalTB(0), nCorrectTBHierarchy(0), nTotalTBHierarchy(0), nCorrectCR(0), nTotalCR(0);
    int nFakeTBHierarchy(0), nFakeCR(0), nSplitTBHierarchy(0), nSplitCR(0), nLost(0), mcPrimaryIndex(0), nTargetMatches(0), nTargetTBHierarchyMatches(0);
    int nTargetCRMatches(0), nTargetGoodTBHierarchyMatches(0), nTargetTBHierarchySplits(0), nTargetTBHierarchyLosses(0);
    IntVector mcPrimaryId, mcPrimaryPdg, mcPrimaryTier, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedTBHierarchyPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoTier, bestMatchPfoIsTestBeam;
    IntVector bestMatchPfoRecoTBId, bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;

    std::stringstream targetSS;
    const std::string name("TB");

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!isTargetPrimary)
            continue;

        // Parent in hierarchy needed even if no match
        const bool hasVisibleTargets((!triggeredToLeading.empty() && LArMCParticleHelper::IsBeamParticle(LArMCParticleHelper::GetParentMCParticle(pMCPrimary))) ?
            triggeredToLeading.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)) != 1 : false);

        if (!hasMatch && !hasVisibleTargets)
            continue;

        associatedMCPrimaries.push_back(pMCPrimary);
        const int nTargetPrimaries(associatedMCPrimaries.size());
        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        const int isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));

        // Leading beam particle is the primary beam particle or a daughter of that particle
        const int isLeadingBeamParticle(LArMCParticleHelper::IsLeadingBeamParticle(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));

        // Tier (0) : Primary, (1) : Daughter, (2) : Granddaughter etc...  Note tier increases for both visible and invisible particles
        const int mcHierarchyTier(LArMCParticleHelper::GetHierarchyTier(pMCPrimary));

        // Identify the number of matched leading particles and flag whether last particle in hierarchy is being considered
        bool isLastTestBeamLeading(false);
        if (isLeadingBeamParticle)
        {
            triggeredToLeadingCounter.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary))++;
            const int nHierarchyLeading(triggeredToLeadingCounter.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
            isLastTestBeamLeading = (nHierarchyLeading == triggeredToLeading.at(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        }

#ifdef MONITORING
        const CartesianVector &targetVertex(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetVertex());
        const float targetVertexX(targetVertex.GetX()), targetVertexY(targetVertex.GetY()), targetVertexZ(targetVertex.GetZ());
#endif

        for (int tier = 0; tier < mcHierarchyTier; tier++) targetSS << " -> ";

        targetSS << (!isTargetPrimary ? "(Non target) " : "")
                 << "PrimaryId " << mcPrimaryIndex
                 << ", TB " << isBeamParticle
                 << ", TB Hierarchy " << isLeadingBeamParticle
                 << ", CR " << isCosmicRay
                 << ", MCPDG " << pMCPrimary->GetParticleId()
                 << ", Tier " << mcHierarchyTier
                 << ", Energy " << pMCPrimary->GetEnergy()
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

        int matchIndex(0), nPrimaryMatches(0), nPrimaryTBHierarchyMatches(0), nPrimaryCRMatches(0), nPrimaryGoodTBHierarchyMatches(0), nPrimaryTBHierarchySplits(0);
#ifdef MONITORING
        float recoVertexX(std::numeric_limits<float>::max()), recoVertexY(std::numeric_limits<float>::max()), recoVertexZ(std::numeric_limits<float>::max());
#endif
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoTestBeam(LArPfoHelper::IsTestBeam(pfoToSharedHits.first));
            const bool isRecoTestBeamHierarchy(LArPfoHelper::IsTestBeam(LArPfoHelper::GetParentPfo(pfoToSharedHits.first)));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            // Tier (0) : Primary, (1) : Daughter, (2) : Granddaughter etc...  Note that the tier only increases for visible particle
            const int pfoHierarchyTier(LArPfoHelper::GetHierarchyTier(pfoToSharedHits.first));
            const int pfoId(pfoToIdMap.at(pfoToSharedHits.first));
            const int recoTBId(isRecoTestBeam || isRecoTestBeamHierarchy ? testBeamPfoToIdMap.at(LArPfoHelper::GetParentPfo(pfoToSharedHits.first)) : -1);

            if (0 == matchIndex++)
            {
                bestMatchPfoId.push_back(pfoId);
                bestMatchPfoPdg.push_back(pfoToSharedHits.first->GetParticleId());
                bestMatchPfoTier.push_back(pfoHierarchyTier);
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
                    const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(pfoToSharedHits.first));
                    recoVertexX = pRecoVertex->GetPosition().GetX();
                    recoVertexY = pRecoVertex->GetPosition().GetY();
                    recoVertexZ = pRecoVertex->GetPosition().GetZ();
                }
                catch (const StatusCodeException &) {}
#endif
            }

            if (isGoodMatch) ++nPrimaryMatches;

            // ATTN: In hierarchy mode let TBHierarchyMatches become effective TBHierarchyMatches and treat the same
            if (isRecoTestBeamHierarchy && isGoodMatch) ++nPrimaryTBHierarchyMatches;
            if (!isRecoTestBeamHierarchy && isGoodMatch) ++nPrimaryCRMatches;

            // Account for splitting of test beam particle into separate reconstructed primary pfos
            const Pfo *const pRecoTB(LArPfoHelper::GetParentPfo(pfoToSharedHits.first));
            const bool isSplitRecoTBHierarchy(!recoTestBeamHierarchies.empty() && !recoTestBeamHierarchies.count(pRecoTB));
            if (!isSplitRecoTBHierarchy && isGoodMatch) ++nPrimaryGoodTBHierarchyMatches;
            if (isSplitRecoTBHierarchy && isLeadingBeamParticle && isGoodMatch) ++nPrimaryTBHierarchySplits;
            recoTestBeamHierarchies.insert(pRecoTB);

            for (int tier = 0; tier < mcHierarchyTier; tier++) targetSS << "    ";

            targetSS << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                     << "MatchedPfoId " << pfoId
                     << ", TB " << isRecoTestBeam
                     << ", TB Hierarchy " << isRecoTestBeamHierarchy;
            if (isRecoTestBeamHierarchy) targetSS << " [TBId: " << recoTBId << "]";
            targetSS << ", CR " << (!isRecoTestBeam && !isRecoTestBeamHierarchy)
                     << ", PDG " << pfoToSharedHits.first->GetParticleId()
                     << ", Tier " << pfoHierarchyTier
                     << ", nMatchedHits " << sharedHitList.size()
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
            for (int tier = 0; tier < mcHierarchyTier; tier++) targetSS << "    ";
            targetSS << "-No matched Pfo" << std::endl;
            bestMatchPfoId.push_back(-1); bestMatchPfoPdg.push_back(0); bestMatchPfoTier.push_back(-1);
            bestMatchPfoIsTestBeam.push_back(0); bestMatchPfoRecoTBId.push_back(-1);
            bestMatchPfoNHitsTotal.push_back(0); bestMatchPfoNHitsU.push_back(0); bestMatchPfoNHitsV.push_back(0); bestMatchPfoNHitsW.push_back(0);
            bestMatchPfoNSharedHitsTotal.push_back(0); bestMatchPfoNSharedHitsU.push_back(0); bestMatchPfoNSharedHitsV.push_back(0); bestMatchPfoNSharedHitsW.push_back(0);
        }

        nPrimaryMatchedPfos.push_back(nPrimaryMatches);
        nPrimaryMatchedTBHierarchyPfos.push_back(nPrimaryTBHierarchyMatches);
        nPrimaryMatchedCRPfos.push_back(nPrimaryCRMatches);
        nTargetMatches += nPrimaryMatches;
        nTargetTBHierarchyMatches += nPrimaryTBHierarchyMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodTBHierarchyMatches += nPrimaryGoodTBHierarchyMatches;
        nTargetTBHierarchySplits += nPrimaryTBHierarchySplits;
        if (0 == nPrimaryMatches) ++nTargetTBHierarchyLosses;

	if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedTBHierarchyPfos", &nPrimaryMatchedTBHierarchyPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedCRPfos", &nPrimaryMatchedCRPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoId", &bestMatchPfoId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPdg", &bestMatchPfoPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTier", &bestMatchPfoTier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsTotal", &bestMatchPfoNHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsU", &bestMatchPfoNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsV", &bestMatchPfoNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsW", &bestMatchPfoNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsTotal", &bestMatchPfoNSharedHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsU", &bestMatchPfoNSharedHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsV", &bestMatchPfoNSharedHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsW", &bestMatchPfoNSharedHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetMatches", nTargetMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetTBHierarchyMatches", nTargetTBHierarchyMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetCRMatches", nTargetCRMatches));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsTestBeam", &bestMatchPfoIsTestBeam));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoRecoTBId", &bestMatchPfoRecoTBId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetGoodTBHierarchyMatches", nTargetGoodTBHierarchyMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetTBHierarchySplits", nTargetTBHierarchySplits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetTBHierarchyLosses", nTargetTBHierarchyLosses));
        }

        if (isCosmicRay || isLastTestBeamLeading)
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries, true));
#ifdef MONITORING
            const int interactionTypeInt(static_cast<int>(interactionType));
#endif
            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectTB(isBeamParticle && (nTargetTBHierarchyMatches == 1) && (nTargetCRMatches == 0));
            const int isCorrectTBHierarchy(isLeadingBeamParticle && (nTargetGoodTBHierarchyMatches == nTargetTBHierarchyMatches) && (nTargetGoodTBHierarchyMatches == nTargetPrimaries) && (nTargetCRMatches == 0) && (nTargetTBHierarchySplits == 0) && (nTargetTBHierarchyLosses == 0));
            const int isCorrectCR(isCosmicRay && (nTargetTBHierarchyMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeTBHierarchy(isCosmicRay && (nTargetTBHierarchyMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitTBHierarchy(!isCosmicRay && ((nTargetTBHierarchyMatches > nTargetPrimaries) || (nTargetTBHierarchySplits > 0)));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            std::stringstream outcomeSS;
            const bool isBeamHierarchy((mcNuanceCode == 2001) | (mcNuanceCode == 2000));
            outcomeSS << LArInteractionTypeHelper::ToString(interactionType) << " (Nuance " << mcNuanceCode << ", TB " << isBeamHierarchy << ", CR " << isCosmicRay << ")" << std::endl;

            if (isBeamParticle) ++nTotalTB;
            if (isLastTestBeamLeading) ++nTotalTBHierarchy;
            if (isCosmicRay) ++nTotalCR;
            if (isCorrectTB) ++nCorrectTB;
            if (isCorrectTBHierarchy) ++nCorrectTBHierarchy;
            if (isCorrectCR) ++nCorrectCR;
            if (isFakeTBHierarchy) ++nFakeTBHierarchy;
            if (isFakeCR) ++nFakeCR;
            if (isSplitTBHierarchy) ++nSplitTBHierarchy;
            if (isSplitCR) ++nSplitCR;
            if (isLost) ++nLost;

            if (isCorrectTBHierarchy) outcomeSS << "IsCorrectTBHierarchy";
            if (isCorrectCR) outcomeSS << "IsCorrectCR ";
            if (isFakeTBHierarchy) outcomeSS << "IsFake" << name << " ";
            if (isFakeCR) outcomeSS << "IsFakeCR ";
            if (isSplitTBHierarchy) outcomeSS << "isSplit" << name << " ";
            if (isSplitCR) outcomeSS << "IsSplitCR ";
            if (isLost) outcomeSS << "IsLost ";
            if (nTargetTBHierarchyMatches > 0) outcomeSS << "(N" << name << "Matches: " << nTargetTBHierarchyMatches << ") ";
            if (nTargetTBHierarchyLosses > 0) outcomeSS << "(N" << name << "Losses: " << nTargetTBHierarchyLosses << ") ";
            if (nTargetTBHierarchySplits > 0) outcomeSS << "(N" << name << "Splits: " << nTargetTBHierarchySplits << ") ";
            if (nTargetCRMatches > 0) outcomeSS << "(NCRMatches: " << nTargetCRMatches << ") ";
            if (printToScreen) std::cout << outcomeSS.str() << std::endl << targetSS.str() << std::endl;

            if (fillTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionTypeInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectTBHierarchy", isCorrectTBHierarchy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectCR", isCorrectCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeTBHierarchy", isFakeTBHierarchy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeCR", isFakeCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitTBHierarchy", isSplitTBHierarchy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitCR", isSplitCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isLost", isLost));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }

            targetSS.str(std::string()); targetSS.clear();
            associatedMCPrimaries.clear();
            nTargetMatches = 0; nTargetTBHierarchyMatches = 0; nTargetCRMatches = 0; nTargetGoodTBHierarchyMatches = 0; nTargetTBHierarchySplits = 0; nTargetTBHierarchyLosses = 0;
            mcPrimaryId.clear(); mcPrimaryPdg.clear(); mcPrimaryTier.clear(); nMCHitsTotal.clear(); nMCHitsU.clear(); nMCHitsV.clear(); nMCHitsW.clear();
            mcPrimaryE.clear(); mcPrimaryPX.clear(); mcPrimaryPY.clear(); mcPrimaryPZ.clear();
            mcPrimaryVtxX.clear(); mcPrimaryVtxY.clear(); mcPrimaryVtxZ.clear(); mcPrimaryEndX.clear(); mcPrimaryEndY.clear(); mcPrimaryEndZ.clear();
            nPrimaryMatchedPfos.clear(); nPrimaryMatchedTBHierarchyPfos.clear(); nPrimaryMatchedCRPfos.clear();
            bestMatchPfoId.clear(); bestMatchPfoPdg.clear(); bestMatchPfoTier.clear(); bestMatchPfoIsTestBeam.clear(); bestMatchPfoRecoTBId.clear();
            bestMatchPfoNHitsTotal.clear(); bestMatchPfoNHitsU.clear(); bestMatchPfoNHitsV.clear(); bestMatchPfoNHitsW.clear();
            bestMatchPfoNSharedHitsTotal.clear(); bestMatchPfoNSharedHitsU.clear(); bestMatchPfoNSharedHitsV.clear(); bestMatchPfoNSharedHitsW.clear();
        }
    }

    if (useInterpretedMatching)
    {
        std::stringstream summarySS;
        summarySS << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalTBHierarchy > 0) summarySS << "#CorrectTBHierarchy: " << nCorrectTBHierarchy << "/" << nTotalTBHierarchy << ", Fraction: " << (nTotalTBHierarchy > 0 ? static_cast<float>(nCorrectTBHierarchy) / static_cast<float>(nTotalTBHierarchy) : 0.f) << std::endl;
        if (nTotalCR > 0) summarySS << "#CorrectCR: " << nCorrectCR << "/" << nTotalCR << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if (nFakeTBHierarchy > 0) summarySS << "#Fake" << name << ": " << nFakeTBHierarchy << " ";
        if (nFakeCR > 0) summarySS << "#FakeCR: " << nFakeCR << " ";
        if (nSplitTBHierarchy > 0) summarySS << "#Split" << name << ": " << nSplitTBHierarchy << " ";
        if (nSplitCR > 0) summarySS << "#SplitCR: " << nSplitCR << " ";
        if (nLost > 0) summarySS << "#Lost: " << nLost << " ";
        if (nFakeTBHierarchy || nFakeCR || nSplitTBHierarchy || nSplitCR || nLost) summarySS << std::endl;
        if (printToScreen) std::cout << summarySS.str();
    }

    if (printToScreen) std::cout << "------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamHierarchyEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
