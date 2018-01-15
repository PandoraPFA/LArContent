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
//    m_printMatchingToScreen(false),
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
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, {mcParticleToHitsMap}, pfoToHitsMap);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {mcParticleToHitsMap}, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Print raw matching information to terminal
    if (m_printAllToScreen)
        this->PrintAllOutput(mcParticleToHitsMap, pfoToHitsMap, mcParticleToPfoHitSharingMap);

//    // Write raw matching information to root file
//    if (m_writeToTree)
//        this->WriteAllOutput(mcParticleToPfoHitSharingMap);
//
//    if (m_printMatchingToScreen || m_visualizeMatching)
//    {
//        // Obtain map: [simple mc primary -> interpreted list of simple matched pfos]
//        MatchingDetailsMap matchingDetailsMap;
//        this->PerformMatching(mcParticleToPfoHitSharingMap, matchingDetailsMap);
//
//        // Print interpreted matching information to terminal
//        if (m_printMatchingToScreen)
//            this->PrintMatchingOutput(mcParticleToPfoHitSharingMap, matchingDetailsMap);
//#ifdef MONITORING
//        // Visualize interpreted matching information
//        if (m_visualizeMatching)
//            this->VisualizeMatchingOutput(mcParticleToPfoHitSharingMap, matchingDetailsMap);
//#endif
//    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintAllOutput(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const
{
    std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;
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

        LArMCParticleHelper::PfoToSharedHitsVector pfoToSharedHitsVector(mcParticleToPfoHitSharingMap.at(pMCPrimary));
        std::sort(pfoToSharedHitsVector.begin(), pfoToSharedHitsVector.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool
        {
            if (a.second.size() != b.second.size())
                return (a.second.size() > b.second.size());

            return LArPfoHelper::SortByNHits(a.first, b.first);
        });

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : pfoToSharedHitsVector)
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
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
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
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::PerformMatching(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, MatchingDetailsMap &matchingDetailsMap) const
//{
//    // Get best matches, one-by-one, until no more strong matches possible
//    IntSet usedMCIds, usedPfoIds;
//    while (GetStrongestPfoMatch(mcPrimaryMatchingMap, usedMCIds, usedPfoIds, matchingDetailsMap)) {}
//
//    // Assign any remaining pfos to primaries, based on number of matched hits
//    GetRemainingPfoMatches(mcPrimaryMatchingMap, usedPfoIds, matchingDetailsMap);
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::GetStrongestPfoMatch(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, IntSet &usedMCIds, IntSet &usedPfoIds,
//    MatchingDetailsMap &matchingDetailsMap) const
//{
//    int bestPfoMatchId(-1);
//    MatchingDetails bestMatchingDetails;
//// switch to new constructs, with addresses rather than ids
//    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
//    {
//        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
//
//        // difficult to maintain this functionality - could run with no thresholds, then recheck at this point?
//        if (!m_useSmallPrimaries && !this->IsGoodMCPrimary(simpleMCPrimary))
//            continue;
//
//        if (usedMCIds.count(simpleMCPrimary.m_id))
//            continue;
//
//        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
//        {
//            if (usedPfoIds.count(simpleMatchedPfo.m_id))
//                continue;
//
//            if (!this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo))
//                continue;
//
//            if (simpleMatchedPfo.m_nMatchedHitsTotal > bestMatchingDetails.m_nMatchedHits)
//            {
//                bestPfoMatchId = simpleMatchedPfo.m_id;
//                bestMatchingDetails.m_matchedPrimaryId = simpleMCPrimary.m_id;
//                bestMatchingDetails.m_nMatchedHits = simpleMatchedPfo.m_nMatchedHitsTotal;
//            }
//        }
//    }
//
//    if (bestPfoMatchId > -1)
//    {
//        matchingDetailsMap[bestPfoMatchId] = bestMatchingDetails;
//        usedMCIds.insert(bestMatchingDetails.m_matchedPrimaryId);
//        usedPfoIds.insert(bestPfoMatchId);
//        return true;
//    }
//
//    return false;
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::GetRemainingPfoMatches(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const IntSet &usedPfoIds,
//    MatchingDetailsMap &matchingDetailsMap) const
//{
//    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
//    {
//        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
//
//        if (!m_useSmallPrimaries && !this->IsGoodMCPrimary(simpleMCPrimary))
//            continue;
//
//        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
//        {
//            if (usedPfoIds.count(simpleMatchedPfo.m_id))
//                continue;
//
//            MatchingDetails &matchingDetails(matchingDetailsMap[simpleMatchedPfo.m_id]);
//
//            if (simpleMatchedPfo.m_nMatchedHitsTotal > matchingDetails.m_nMatchedHits)
//            {
//                matchingDetails.m_matchedPrimaryId = simpleMCPrimary.m_id;
//                matchingDetails.m_nMatchedHits = simpleMatchedPfo.m_nMatchedHitsTotal;
//            }
//        }
//    }
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::PrintMatchingOutput(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const
//{
//    std::cout << "---PROCESSED-MATCHING-OUTPUT--------------------------------------------------------------------" << std::endl;
//    std::cout << "MinPrimaryGoodHits " << m_matchingMinPrimaryHits << ", MinHitsForGoodView " << m_matchingMinHitsForGoodView << ", MinPrimaryGoodViews " << m_matchingMinPrimaryGoodViews << std::endl;
//    std::cout << "UseSmallPrimaries " << m_useSmallPrimaries << ", MinSharedHits " << m_matchingMinSharedHits << ", MinCompleteness " << m_matchingMinCompleteness << ", MinPurity " << m_matchingMinPurity << std::endl;
//
//    bool isCorrect(true), isCalculable(false);
//
//    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
//    {
//        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
//        const bool hasMatch(this->HasMatch(simpleMCPrimary, mapValue.second, matchingDetailsMap));
//        const bool isTargetPrimary(this->IsGoodMCPrimary(simpleMCPrimary) && (NEUTRON != simpleMCPrimary.m_pdgCode));
//
//        if (!hasMatch && !isTargetPrimary)
//            continue;
//
//        std::cout << std::endl << (!isTargetPrimary ? "(Non target) " : "")
//                  << "Primary " << simpleMCPrimary.m_id << ", PDG " << simpleMCPrimary.m_pdgCode << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
//                  << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << "),"
//                  << " [nGood " << simpleMCPrimary.m_nGoodMCHitsTotal << " (" << simpleMCPrimary.m_nGoodMCHitsU << ", " << simpleMCPrimary.m_nGoodMCHitsV
//                  << ", " << simpleMCPrimary.m_nGoodMCHitsW << ")]" << std::endl;
//
//        if (NEUTRON != simpleMCPrimary.m_pdgCode)
//            isCalculable = true;
//
//        unsigned int nMatches(0);
//
//        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
//        {
//            if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
//            {
//                const bool isGoodMatch(this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo));
//
//                if (isGoodMatch) ++nMatches;
//                std::cout << "-" << (!isGoodMatch ? "(Below threshold) " : "") << "MatchedPfo " << simpleMatchedPfo.m_id;
//
//                if (simpleMatchedPfo.m_parentId >= 0) std::cout << ", ParentPfo " << simpleMatchedPfo.m_parentId;
//
//                std::cout << ", PDG " << simpleMatchedPfo.m_pdgCode << ", nMatchedHits " << simpleMatchedPfo.m_nMatchedHitsTotal
//                          << " (" << simpleMatchedPfo.m_nMatchedHitsU << ", " << simpleMatchedPfo.m_nMatchedHitsV << ", " << simpleMatchedPfo.m_nMatchedHitsW << ")"
//                          << ", nPfoHits " << simpleMatchedPfo.m_nPfoHitsTotal
//                          << " (" << simpleMatchedPfo.m_nPfoHitsU << ", " << simpleMatchedPfo.m_nPfoHitsV << ", " << simpleMatchedPfo.m_nPfoHitsW << ")" << std::endl;
//            }
//        }
//
//        if (isTargetPrimary && (1 != nMatches))
//            isCorrect = false;
//    }
//
//    // is correct from neutrino point of view, from beam particle view and from pov of each cr...
//    std::cout << std::endl << "Is correct? " << (isCorrect && isCalculable) << std::endl;
//    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//#ifdef MONITORING
//void EventValidationAlgorithm::VisualizeMatchingOutput(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector,
//    const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const
//{
//    std::cout << "---VISUALIZE-MATCHING-OUTPUT--------------------------------------------------------------------" << std::endl;
//    std::cout << "MinPrimaryGoodHits " << m_matchingMinPrimaryHits << ", MinHitsForGoodView " << m_matchingMinHitsForGoodView << ", MinPrimaryGoodViews " << m_matchingMinPrimaryGoodViews << std::endl;
//    std::cout << "UseSmallPrimaries " << m_useSmallPrimaries << ", MinSharedHits " << m_matchingMinSharedHits << ", MinCompleteness " << m_matchingMinCompleteness << ", MinPurity " << m_matchingMinPurity << std::endl;
//
//    // needs to work for neutrinos, beam particles and crs too
//    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), m_visualizeGaps, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
//    const HitTypeVector hitTypeVector{TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};
//
//    for (const HitType hitType : hitTypeVector)
//    {
//        PfoSet allPrimaryMatchedPfos;
//        const std::string hitTypeString((TPC_VIEW_U == hitType) ? "_U" : (TPC_VIEW_V == hitType) ? "_V" : "_W");
//
//        if (m_visualizeVertices)
//            this->VisualizeVertexMatches(mcNeutrinoVector, recoNeutrinoVector, hitType);
//
//        for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
//        {
//            const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
//            bool hasMatch(this->HasMatch(simpleMCPrimary, mapValue.second, matchingDetailsMap));
//            const bool isTargetPrimary(this->IsGoodMCPrimary(simpleMCPrimary) && (NEUTRON != simpleMCPrimary.m_pdgCode));
//
//            if (!hasMatch && !isTargetPrimary)
//                continue;
//
//            PfoSet belowThresholdPfos;
//            PfoVector primaryMatchedPfos;
//
//            for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
//            {
//                if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
//                {
//                    primaryMatchedPfos.push_back(simpleMatchedPfo.m_pPandoraAddress);
//
//                    if (!this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo))
//                        belowThresholdPfos.insert(simpleMatchedPfo.m_pPandoraAddress);
//                }
//            }
//
//            std::string name("UNKNOWN_"); Color color(BLACK);
//            this->GetPrimaryDetails(simpleMCPrimary, mcPrimaryMatchingMap, name, color);
//            name.insert(0, !isTargetPrimary ? "NonTarget" : "");
//            name += hitTypeString;
//
//            if (isTargetPrimary && primaryMatchedPfos.empty())
//            {
//                const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), simpleMCPrimary.m_vertex, hitType));
//                const CartesianVector endpointPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), simpleMCPrimary.m_endpoint, hitType));
//                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertexPosition2D, "MissingPrimaryVtx_" + name, RED, 1);
//                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpointPosition2D, "MissingPrimaryEnd_" + name, RED, 1);
//            }
//
//            for (const ParticleFlowObject *const pPrimaryPfo : primaryMatchedPfos)
//            {
//                const bool isSplit((static_cast<int>(primaryMatchedPfos.size()) - static_cast<int>(belowThresholdPfos.size())) > 1);
//                const bool isBestMatch(pPrimaryPfo == primaryMatchedPfos.front());
//                const bool isBelowThreshold(belowThresholdPfos.count(pPrimaryPfo));
//                const std::string prefix(isBelowThreshold ? "BelowThreshold" : !isSplit ? "Matched" : isBestMatch ? "BestFragment" : "Fragment");
//
//                PfoList allPfos;
//                LArPfoHelper::GetAllDownstreamPfos(pPrimaryPfo, allPfos);
//                PfoVector allSortedPfos(allPfos.begin(), allPfos.end());
//                std::sort(allSortedPfos.begin(), allSortedPfos.end(), LArPfoHelper::SortByNHits);
//
//                for (const ParticleFlowObject *const pPfo : allSortedPfos)
//                {
//                    ClusterList clusterList;
//                    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);
//                    const std::string hierarchy((pPfo != pPrimaryPfo) ? "Daughter_" : "");
//                    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &clusterList, prefix + "Clusters_" + hierarchy + name, color);
//
//                    if (!isBelowThreshold && isSplit && !isBestMatch && !clusterList.empty())
//                    {
//                        CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
//                        LArClusterHelper::GetExtremalCoordinates(*(clusterList.begin()), innerCoordinate, outerCoordinate);
//                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerCoordinate, "SplitMarker_" + hierarchy + name, MAGENTA, 1);
//                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerCoordinate, "SplitMarker_" + hierarchy + name, MAGENTA, 1);
//                    }
//                }
//            }
//
//            allPrimaryMatchedPfos.insert(primaryMatchedPfos.begin(), primaryMatchedPfos.end());
//        }
//
//        if (m_visualizeRemnants)
//            this->VisualizeRemnants(hitType);
//
//        PandoraMonitoringApi::ViewEvent(this->GetPandora());
//    }
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::VisualizeVertexMatches(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector, const HitType hitType) const
//{
//    const Vertex *pBestVertex(nullptr);
//    float closestDistance(m_vertexVisualizationDeltaR);
//
//    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
//    {
//        for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinoVector)
//        {
//            const Vertex *const pVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
//            const float distance((pMCNeutrino->GetEndpoint() - pVertex->GetPosition()).GetMagnitude());
//
//            if (distance < closestDistance)
//            {
//                closestDistance = distance;
//                pBestVertex = pVertex;
//            }
//        }
//    }
//
//    const std::string hitTypeString((TPC_VIEW_U == hitType) ? "_U" : (TPC_VIEW_V == hitType) ? "_V" : "_W");
//    const bool isGood((1 == mcNeutrinoVector.size()) && (1 == recoNeutrinoVector.size()) && pBestVertex);
//    const Color color(isGood ? GRAY : ORANGE);
//
//    unsigned int displayIndex(0);
//
//    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
//    {
//        const CartesianVector mcVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCNeutrino->GetEndpoint(), hitType));
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertex2D, "MCNeutrinoVertex_" + TypeToString(displayIndex++) + hitTypeString, color, 1);
//    }
//
//    displayIndex = 0;
//
//    for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinoVector)
//    {
//        const Vertex *const pThisVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
//        const bool isThisGood(isGood || ((1 == mcNeutrinoVector.size()) && (pThisVertex == pBestVertex)));
//        const std::string thisPrefix(isThisGood ? "Good" : "Displaced");
//        const Color thisColor(isThisGood ? CYAN : DARKVIOLET);
//
//        const CartesianVector recoVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pThisVertex->GetPosition(), hitType));
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &recoVertex2D, thisPrefix + "RecoNeutrinoVertex_" + TypeToString(displayIndex++) + hitTypeString, thisColor, 1);
//    }
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::GetPrimaryDetails(const SimpleMCPrimary &thisSimpleMCPrimary, const MCPrimaryMatchingMap &mcPrimaryMatchingMap,
//    std::string &name, Color &color) const
//{
//    // ATTN: Relies on fact that mcPrimaryMatchingMap is sorted by number of true hits
//    unsigned int nLikeParticles(0);
//
//    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
//    {
//        const SimpleMCPrimary &mapSimpleMCPrimary(mapValue.first);
//
//        if (thisSimpleMCPrimary.m_pdgCode == mapSimpleMCPrimary.m_pdgCode)
//            ++nLikeParticles;
//
//        if (thisSimpleMCPrimary.m_id == mapSimpleMCPrimary.m_id)
//        {
//            const std::string index(TypeToString(nLikeParticles));
//
//            switch (mapSimpleMCPrimary.m_pdgCode)
//            {
//                case MU_MINUS: {name = "MU_MINUS" + index; color = RED;     return;}
//                case MU_PLUS:  {name = "MU_PLUS"  + index; color = RED;     return;}
//                case E_MINUS:  {name = "E_MINUS"  + index; color = ORANGE;  return;}
//                case E_PLUS:   {name = "E_PLUS"   + index; color = ORANGE;  return;}
//                case PROTON:   {name = "PROTON"   + index; color = (nLikeParticles % 2 == 0) ? VIOLET : BLUE;  return;}
//                case PI_MINUS: {name = "PI_MINUS" + index; color = MAGENTA; return;}
//                case PI_PLUS:  {name = "PI_PLUS"  + index; color = MAGENTA; return;}
//                case PHOTON:   {name = "PHOTON"   + index; color = (nLikeParticles % 2 == 0) ? TEAL   : GREEN; return;}
//                case NEUTRON:  {name = "NEUTRON"  + index; color = CYAN;    return;}
//                default: return;
//            }
//        }
//    }
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::VisualizeRemnants(const HitType hitType) const
//{
//    const std::string hitTypeString((TPC_VIEW_U == hitType) ? "_U" : (TPC_VIEW_V == hitType) ? "_V" : "_W");
//
//    const CaloHitList *pCaloHitList = nullptr;
//    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
//
//    CaloHitList availableHits;
//
//    for (const CaloHit *const pCaloHit : *pCaloHitList)
//    {
//        if (PandoraContentApi::IsAvailable(*this, pCaloHit) && (hitType == pCaloHit->GetHitType()))
//            availableHits.push_back(pCaloHit);
//    }
//
//    if (!availableHits.empty())
//        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &availableHits, "AllRemnantHits" + hitTypeString, LIGHTCYAN);
//
//    CaloHitList isolatedHits;
//    ClusterList availableClusters;
//
//    for (const std::string &clusterListName : m_clusterListNames)
//    {
//        const ClusterList *pClusterList = nullptr;
//        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
//
//        if (!pClusterList)
//            continue;
//
//        for (const Cluster *const pCluster : *pClusterList)
//        {
//            if (PandoraContentApi::IsAvailable(*this, pCluster) && (hitType == LArClusterHelper::GetClusterHitType(pCluster)))
//                availableClusters.push_back(pCluster);
//
//            if (!PandoraContentApi::IsAvailable(*this, pCluster) && (hitType == LArClusterHelper::GetClusterHitType(pCluster)))
//                isolatedHits.insert(isolatedHits.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
//        }
//    }
//
//    if (!availableClusters.empty())
//        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &availableClusters, "AllRemnantClusters" + hitTypeString, GRAY);
//
//    if (!isolatedHits.empty())
//        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &isolatedHits, "AllIsolatedHitsInParticles" + hitTypeString, LIGHTYELLOW);
//}
//#endif
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::IsGoodMCPrimary(const SimpleMCPrimary &simpleMCPrimary) const
//{
//    if (simpleMCPrimary.m_nGoodMCHitsTotal < m_matchingMinPrimaryHits)
//        return false;
//
//    int nGoodViews(0);
//    if (simpleMCPrimary.m_nGoodMCHitsU >= m_matchingMinHitsForGoodView) ++nGoodViews;
//    if (simpleMCPrimary.m_nGoodMCHitsV >= m_matchingMinHitsForGoodView) ++nGoodViews;
//    if (simpleMCPrimary.m_nGoodMCHitsW >= m_matchingMinHitsForGoodView) ++nGoodViews;
//
//    if (nGoodViews < m_matchingMinPrimaryGoodViews)
//        return false;
//
//    return true;
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::HasMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfoList &simpleMatchedPfoList,
//    const MatchingDetailsMap &matchingDetailsMap) const
//{
//    for (const SimpleMatchedPfo &simpleMatchedPfo : simpleMatchedPfoList)
//    {
//        if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
//            return true;
//    }
//
//    return false;
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::IsGoodMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfo &simpleMatchedPfo) const
//{
//    const float purity((simpleMatchedPfo.m_nPfoHitsTotal > 0) ? static_cast<float>(simpleMatchedPfo.m_nMatchedHitsTotal) / static_cast<float>(simpleMatchedPfo.m_nPfoHitsTotal) : 0.f);
//    const float completeness((simpleMCPrimary.m_nMCHitsTotal > 0) ? static_cast<float>(simpleMatchedPfo.m_nMatchedHitsTotal) / static_cast<float>(simpleMCPrimary.m_nMCHitsTotal) : 0.f);
//
//    return ((simpleMatchedPfo.m_nMatchedHitsTotal >= m_matchingMinSharedHits) && (purity >= m_matchingMinPurity) && (completeness >= m_matchingMinCompleteness));
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//void EventValidationAlgorithm::GetPfoIdMap(const PfoList &pfoList, PfoIdMap &pfoIdMap) const
//{
//    PfoVector pfoVector(pfoList.begin(), pfoList.end());
//    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
//
//    int id(0);
//
//    for (const ParticleFlowObject *const pPfo : pfoVector)
//        (void) pfoIdMap.insert(PfoIdMap::value_type(pPfo, id++));
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::SortSimpleMCPrimaries(const SimpleMCPrimary &lhs, const SimpleMCPrimary &rhs)
//{
//    if (lhs.m_nGoodMCHitsTotal != rhs.m_nGoodMCHitsTotal)
//        return (lhs.m_nGoodMCHitsTotal > rhs.m_nGoodMCHitsTotal);
//
//    return (lhs.m_energy > rhs.m_energy);
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::SortSimpleMatchedPfos(const SimpleMatchedPfo &lhs, const SimpleMatchedPfo &rhs)
//{
//    if (lhs.m_nMatchedHitsTotal != rhs.m_nMatchedHitsTotal)
//        return (lhs.m_nMatchedHitsTotal > rhs.m_nMatchedHitsTotal);
//
//    if (lhs.m_nPfoHitsTotal != rhs.m_nPfoHitsTotal)
//        return (lhs.m_nPfoHitsTotal > rhs.m_nPfoHitsTotal);
//
//    return (lhs.m_id < rhs.m_id);
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::SortRecoNeutrinos(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
//{
//    if (!LArPfoHelper::IsNeutrino(pLhs) || !LArPfoHelper::IsNeutrino(pRhs))
//        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
//
//    PfoList downstreamPfosLhs, downstreamPfosRhs;
//    LArPfoHelper::GetAllDownstreamPfos(pLhs, downstreamPfosLhs);
//    LArPfoHelper::GetAllDownstreamPfos(pRhs, downstreamPfosRhs);
//
//    // If just left with the neutrino pfos themselves
//    if ((1 == downstreamPfosLhs.size()) && (1 == downstreamPfosRhs.size()))
//    {
//        // ATTN Not a good pair of tie-breakers, but this should be rare (ideally shouldn't have any neutrinos without daughter pfos)
//        if (!pLhs->GetVertexList().empty() && !pRhs->GetVertexList().empty())
//            return ((*(pLhs->GetVertexList().begin()))->GetPosition().GetZ() < (*(pRhs->GetVertexList().begin()))->GetPosition().GetZ());
//
//        return (pLhs->GetParticleId() < pRhs->GetParticleId());
//    }
//
//    PfoVector pfoVectorLhs(downstreamPfosLhs.begin(), downstreamPfosLhs.end());
//    PfoVector pfoVectorRhs(downstreamPfosRhs.begin(), downstreamPfosRhs.end());
//    std::sort(pfoVectorLhs.begin(), pfoVectorLhs.end(), LArPfoHelper::SortByNHits);
//    std::sort(pfoVectorRhs.begin(), pfoVectorRhs.end(), LArPfoHelper::SortByNHits);
//
//    if (pfoVectorLhs.empty() || pfoVectorRhs.empty())
//        throw StatusCodeException(STATUS_CODE_FAILURE);
//
//    return (LArPfoHelper::SortByNHits(pfoVectorLhs.front(), pfoVectorRhs.front()));
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------------------------------------
//
//EventValidationAlgorithm::SimpleMCPrimary::SimpleMCPrimary() :
//    m_id(-1),
//    m_pdgCode(0),
//    m_nMCHitsTotal(0),
//    m_nMCHitsU(0),
//    m_nMCHitsV(0),
//    m_nMCHitsW(0),
//    m_nGoodMCHitsTotal(0),
//    m_nGoodMCHitsU(0),
//    m_nGoodMCHitsV(0),
//    m_nGoodMCHitsW(0),
//    m_energy(0.f),
//    m_momentum(0.f, 0.f, 0.f),
//    m_vertex(-1.f, -1.f, -1.f),
//    m_endpoint(-1.f, -1.f, -1.f),
//    m_nMatchedPfos(0),
//    m_pPandoraAddress(nullptr)
//{
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
//
//bool EventValidationAlgorithm::SimpleMCPrimary::operator<(const SimpleMCPrimary &rhs) const
//{
//    if (this == &rhs)
//        return false;
//
//    return (m_id < rhs.m_id);
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------------------------------------
//
//EventValidationAlgorithm::SimpleMatchedPfo::SimpleMatchedPfo() :
//    m_id(-1),
//    m_parentId(-1),
//    m_pdgCode(0),
//    m_nPfoHitsTotal(0),
//    m_nPfoHitsU(0),
//    m_nPfoHitsV(0),
//    m_nPfoHitsW(0),
//    m_nMatchedHitsTotal(0),
//    m_nMatchedHitsU(0),
//    m_nMatchedHitsV(0),
//    m_nMatchedHitsW(0),
//    m_vertex(0.f, 0.f, 0.f),
//    m_endpoint(0.f, 0.f, 0.f),
//    m_vertexDirection(0.f, 0.f, 0.f),
//    m_endDirection(0.f, 0.f, 0.f),
//    m_pPandoraAddress(nullptr)
//{
//}
//
////------------------------------------------------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------------------------------------------------
//
//EventValidationAlgorithm::MatchingDetails::MatchingDetails() :
//    m_matchedPrimaryId(-1),
//    m_nMatchedHits(0)
//{
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

//    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
//        "PrintMatchingToScreen", m_printMatchingToScreen));
//
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
