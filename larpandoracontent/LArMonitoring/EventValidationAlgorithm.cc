/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

#include <iomanip>

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
//    m_integrateOverRecoNeutrinos(false),
//    m_useRecoNeutrinosOnly(true),
    m_useTrueNeutrinosOnly(false),
//    m_primaryPfosOnly(true),
//    m_collapseToPrimaryPfos(true),
//    m_selectInputHits(true),
//    m_minHitNeutrinoWeight(0.1f),
//    m_minHitSharingFraction(0.9f),
//    m_maxPhotonPropagation(2.5f),
    m_printAllToScreen(false),
    m_printMatchingToScreen(true),
//    m_visualizeMatching(false),
//    m_visualizeVertices(false),
//    m_visualizeRemnants(false),
//    m_visualizeGaps(false),
    m_writeToTree(false),
//    m_matchingMinPrimaryHits(15),
//    m_matchingMinHitsForGoodView(5),
//    m_matchingMinPrimaryGoodViews(2),
//    m_useSmallPrimaries(true),
//    m_matchingMinSharedHits(5),
//    m_matchingMinCompleteness(0.1f),
//    m_matchingMinPurity(0.5f),
//    m_vertexVisualizationDeltaR(1.f),
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
    parameters.m_minHitSharingFraction = 0.f;

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
        this->PrintOutput(mcParticleToHitsMap, pfoToHitsMap, mcParticleToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMCToPfoHitSharingMap(mcParticleToHitsMap, mcParticleToPfoHitSharingMap, interpretedMCToPfoHitSharingMap);

    if (m_printMatchingToScreen)
        this->PrintOutput(mcParticleToHitsMap, pfoToHitsMap, interpretedMCToPfoHitSharingMap);

//    if (m_writeToTree)
//        this->WriteOutput(mcParticleToHitsMap, pfoToHitsMap, interpretedMCToPfoHitSharingMap);
//#ifdef MONITORING
//    if (m_visualizeMatching)
//        this->VisualizeOutput(mcParticleToHitsMap, pfoToHitsMap, interpretedMCToPfoHitSharingMap);
//#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintOutput(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const
{
    std::cout << "---MATCHING-OUTPUT------------------------------------------------------------------------------" << std::endl;
    // TODO print true neutrinos, vertex info

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcParticleToHitsMap}, mcPrimaryVector);

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToHitsMap, primaryPfoVector);

    PfoToIdMap pfoToIdMap;
    unsigned int primaryIndex(0), pfoIndex(0);

    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
        (void) pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, pfoIndex++));

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const CaloHitList &mcPrimaryHitList(mcParticleToHitsMap.at(pMCPrimary));

        std::cout << std::endl << "--Primary " << primaryIndex++
                  << ", Nu " << LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary)
                  << ", TB " << LArMCParticleHelper::IsBeamParticle(pMCPrimary)
                  << ", CR " << LArMCParticleHelper::IsCosmicRay(pMCPrimary)
                  << ", MCPDG " << pMCPrimary->GetParticleId()
                  << ", Energy " << pMCPrimary->GetEnergy()
                  << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                  << ", nMCHits " << mcPrimaryHitList.size()
                  << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList)
                  << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList)
                  << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcParticleToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(pfoToHitsMap.at(pfoToSharedHits.first));

            std::cout << "-MatchedPfo " << pfoToIdMap.at(pfoToSharedHits.first)
                << ", PDG " << pfoToSharedHits.first->GetParticleId()
                << ", Nu " << LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first)
                << ", CR " << !LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first)
                << ", nMatchedHits " << sharedHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                << ", nPfoHits " << pfoHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;
        }
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::InterpretMCToPfoHitSharingMap(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcParticleToHitsMap}, mcPrimaryVector);

    PfoSet usedPfos;
    while (this->GetStrongestPfoMatch(mcPrimaryVector, mcToPfoHitSharingMap, usedPfos, interpretedMCToPfoHitSharingMap)) {}
    this->GetRemainingPfoMatches(mcPrimaryVector, mcToPfoHitSharingMap, usedPfos, interpretedMCToPfoHitSharingMap);

    // Ensure all primaries have an entry, and sorting is as desired
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        LArMCParticleHelper::PfoToSharedHitsVector &pfoHitPairs(interpretedMCToPfoHitSharingMap[pMCPrimary]);
        std::sort(pfoHitPairs.begin(), pfoHitPairs.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool {
            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first)); });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::GetStrongestPfoMatch(const MCParticleVector &mcPrimaryVector, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    const MCParticle *pBestMCParticle(nullptr);
    LArMCParticleHelper::PfoCaloHitListPair bestPfoHitPair(nullptr, CaloHitList());

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
//        if (!m_useSmallPrimaries && !this->IsGoodMCPrimary(simpleMCPrimary)) TODO
//            continue;

        if (interpretedMCToPfoHitSharingMap.count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

//            if (!this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo)) TODO
//                continue;

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

void EventValidationAlgorithm::GetRemainingPfoMatches(const MCParticleVector &mcPrimaryVector, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
//        if (!m_useSmallPrimaries && !this->IsGoodMCPrimary(simpleMCPrimary)) TODO
//            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
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

//void EventValidationAlgorithm::WriteAllOutput(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector,
//    const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap, const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const
//{
//#ifdef MONITORING
//    int mcNeutrinoNuance(-1), mcNeutrinoPdg(0);
//    float mcNeutrinoVtxX(-1.f), mcNeutrinoVtxY(-1.f), mcNeutrinoVtxZ(-1.f);
//    float mcNeutrinoE(0.f), mcNeutrinoPX(0.f), mcNeutrinoPY(0.f), mcNeutrinoPZ(0.f);
//    const int nMCNeutrinos(mcNeutrinoVector.size()), nRecoNeutrinos(recoNeutrinoVector.size()), nMCPrimaries(mcPrimaryMatchingMap.size());
//// just get first(?) true neutrino with a nuance code; ought to have n true neutrinos, n true beam particles, n true crs
//    if (!mcNeutrinoVector.empty())
//    {
//        const MCParticle *const pMCNeutrino = mcNeutrinoVector.front();
//
//        mcNeutrinoPdg = pMCNeutrino->GetParticleId();
//        mcNeutrinoVtxX = pMCNeutrino->GetEndpoint().GetX();
//        mcNeutrinoVtxY = pMCNeutrino->GetEndpoint().GetY();
//        mcNeutrinoVtxZ = pMCNeutrino->GetEndpoint().GetZ();
//
//        mcNeutrinoE = pMCNeutrino->GetEnergy();
//        mcNeutrinoPX = pMCNeutrino->GetMomentum().GetX();
//        mcNeutrinoPY = pMCNeutrino->GetMomentum().GetY();
//        mcNeutrinoPZ = pMCNeutrino->GetMomentum().GetZ();
//
//        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);
//
//        if (pLArMCNeutrino)
//            mcNeutrinoNuance = pLArMCNeutrino->GetNuanceCode();
//    }
//
//// TODO reco neutrino or beam particle vertices
//
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCNeutrinos", nMCNeutrinos));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPdg", mcNeutrinoPdg));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoNuance", mcNeutrinoNuance));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nRecoNeutrinos", nRecoNeutrinos));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoPdg", recoNeutrinoPdg));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoVtxX", mcNeutrinoVtxX));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoVtxY", mcNeutrinoVtxY));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoVtxZ", mcNeutrinoVtxZ));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoE", mcNeutrinoE));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPX", mcNeutrinoPX));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPY", mcNeutrinoPY));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPZ", mcNeutrinoPZ));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxX", recoNeutrinoVtxX));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxY", recoNeutrinoVtxY));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxZ", recoNeutrinoVtxZ));
//
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCPrimaries", nMCPrimaries));
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoNPrimaries", recoNeutrinoNPrimaries));
//// TODO n beam particles, n nu primaries, n crs
//
//// now use new structures
//    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
//    {
//        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
//
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryId", simpleMCPrimary.m_id));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPdg", simpleMCPrimary.m_pdgCode));
//// origin: nu, beam ptcle or cr
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsTotal", simpleMCPrimary.m_nMCHitsTotal));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsU", simpleMCPrimary.m_nMCHitsU));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsV", simpleMCPrimary.m_nMCHitsV));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsW", simpleMCPrimary.m_nMCHitsW));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsTotal", simpleMCPrimary.m_nGoodMCHitsTotal));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsU", simpleMCPrimary.m_nGoodMCHitsU));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsV", simpleMCPrimary.m_nGoodMCHitsV));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsW", simpleMCPrimary.m_nGoodMCHitsW));
//
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryE", simpleMCPrimary.m_energy));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPX", simpleMCPrimary.m_momentum.GetX()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPY", simpleMCPrimary.m_momentum.GetY()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPZ", simpleMCPrimary.m_momentum.GetZ()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxX", simpleMCPrimary.m_vertex.GetX()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxY", simpleMCPrimary.m_vertex.GetY()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxZ", simpleMCPrimary.m_vertex.GetZ()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndX", simpleMCPrimary.m_endpoint.GetX()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndY", simpleMCPrimary.m_endpoint.GetY()));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndZ", simpleMCPrimary.m_endpoint.GetZ()));
//
//        const int mcPrimaryNMatchedPfos(mapValue.second.size());
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNMatchedPfos", mcPrimaryNMatchedPfos));
//
//        IntVector pfoIdVector, pfoParentIdVector, pfoPdgVector, pfoNHitsTotalVector, pfoNHitsUVector, pfoNHitsVVector, pfoNHitsWVector,
//            pfoNMatchedHitsTotalVector, pfoNMatchedHitsUVector, pfoNMatchedHitsVVector, pfoNMatchedHitsWVector;
//        FloatVector pfoVtxXVector, pfoVtxYVector, pfoVtxZVector, pfoEndXVector, pfoEndYVector, pfoEndZVector,
//            pfoVtxDirXVector, pfoVtxDirYVector, pfoVtxDirZVector, pfoEndDirXVector, pfoEndDirYVector, pfoEndDirZVector;
//
//        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
//        {
//            pfoIdVector.push_back(simpleMatchedPfo.m_id);
//            pfoParentIdVector.push_back(simpleMatchedPfo.m_parentId);
//            pfoPdgVector.push_back(simpleMatchedPfo.m_pdgCode);
//            pfoNHitsTotalVector.push_back(simpleMatchedPfo.m_nPfoHitsTotal);
//            pfoNHitsUVector.push_back(simpleMatchedPfo.m_nPfoHitsU);
//            pfoNHitsVVector.push_back(simpleMatchedPfo.m_nPfoHitsV);
//            pfoNHitsWVector.push_back(simpleMatchedPfo.m_nPfoHitsW);
//            pfoNMatchedHitsTotalVector.push_back(simpleMatchedPfo.m_nMatchedHitsTotal);
//            pfoNMatchedHitsUVector.push_back(simpleMatchedPfo.m_nMatchedHitsU);
//            pfoNMatchedHitsVVector.push_back(simpleMatchedPfo.m_nMatchedHitsV);
//            pfoNMatchedHitsWVector.push_back(simpleMatchedPfo.m_nMatchedHitsW);
//
//            pfoVtxXVector.push_back(simpleMatchedPfo.m_vertex.GetX());
//            pfoVtxYVector.push_back(simpleMatchedPfo.m_vertex.GetY());
//            pfoVtxZVector.push_back(simpleMatchedPfo.m_vertex.GetZ());
//            pfoEndXVector.push_back(simpleMatchedPfo.m_endpoint.GetX());
//            pfoEndYVector.push_back(simpleMatchedPfo.m_endpoint.GetY());
//            pfoEndZVector.push_back(simpleMatchedPfo.m_endpoint.GetZ());
//            pfoVtxDirXVector.push_back(simpleMatchedPfo.m_vertexDirection.GetX());
//            pfoVtxDirYVector.push_back(simpleMatchedPfo.m_vertexDirection.GetY());
//            pfoVtxDirZVector.push_back(simpleMatchedPfo.m_vertexDirection.GetZ());
//            pfoEndDirXVector.push_back(simpleMatchedPfo.m_endDirection.GetX());
//            pfoEndDirYVector.push_back(simpleMatchedPfo.m_endDirection.GetY());
//            pfoEndDirZVector.push_back(simpleMatchedPfo.m_endDirection.GetZ());
//        }
//// TODO origin reco nu x, beam ptcl y or cr z
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoId", &pfoIdVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoParentId", &pfoParentIdVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoPdg", &pfoPdgVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsTotal", &pfoNHitsTotalVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsU", &pfoNHitsUVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsV", &pfoNHitsVVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsW", &pfoNHitsWVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsTotal", &pfoNMatchedHitsTotalVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsU", &pfoNMatchedHitsUVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsV", &pfoNMatchedHitsVVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsW", &pfoNMatchedHitsWVector));
//
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxX", &pfoVtxXVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxY", &pfoVtxYVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxZ", &pfoVtxZVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndX", &pfoEndXVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndY", &pfoEndYVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndZ", &pfoEndZVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxDirX", &pfoVtxDirXVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxDirY", &pfoVtxDirYVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxDirZ", &pfoVtxDirZVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndDirX", &pfoEndDirXVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndDirY", &pfoEndDirYVector));
//        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndDirZ", &pfoEndDirZVector));
//
//        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
//    }
//#else
//    std::cout << "Monitoring functionality unavailable. nMCNeutrinos " << mcNeutrinoVector.size() << ", nRecoNeutrinos" << recoNeutrinoVector.size()
//              << ", nMCPrimaries " << mcPrimaryMatchingMap.size() << ", nMCParticles " << mcToPrimaryMCMap.size() << std::endl;
//#endif
//}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
         "ClusterListNames", m_clusterListNames));

//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "IntegrateOverRecoNeutrinos", m_integrateOverRecoNeutrinos));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "UseRecoNeutrinosOnly", m_useRecoNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "PrimaryPfosOnly", m_primaryPfosOnly));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "CollapseToPrimaryPfos", m_collapseToPrimaryPfos));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "SelectInputHits", m_selectInputHits));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MinHitNeutrinoWeight", m_minHitNeutrinoWeight));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MinHitSharingFraction", m_minHitSharingFraction));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "VisualizeMatching", m_visualizeMatching));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "VisualizeVertices", m_visualizeVertices));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//         "VisualizeRemnants", m_visualizeRemnants));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "VisualizeGaps", m_visualizeGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MatchingMinPrimaryHits", m_matchingMinPrimaryHits));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MatchingMinHitsForGoodView", m_matchingMinHitsForGoodView));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MatchingMinPrimaryGoodViews", m_matchingMinPrimaryGoodViews));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "UseSmallPrimaries", m_useSmallPrimaries));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MatchingMinSharedHits", m_matchingMinSharedHits));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MatchingMinCompleteness", m_matchingMinCompleteness));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "MatchingMinPurity", m_matchingMinPurity));
//
//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "VertexVisualizationDeltaR", m_vertexVisualizationDeltaR));

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
