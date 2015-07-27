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
#include "LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_printToScreen(true),
    m_writeToTree(true),
    m_minHitsToPrintPrimary(0),
    m_fileIdentifier(0),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::~EventValidationAlgorithm()
{
    if (m_writeToTree)
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
    LArMonitoringHelper::ExtractNeutrinoDaughters(pfoList);

    // Extract monitoring information
    PfoIdMap pfoIdMap;                                              // pfo -> unique identifier
    this->GetPfoIdMap(pfoList, pfoIdMap);

    MCParticleVector mcNeutrinoList;                                // true neutrinos
    LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoList);

    PfoList recoNeutrinoList;                                       // reco neutrinos
    LArMonitoringHelper::GetRecoNeutrinos(pPfoList, recoNeutrinoList);

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

    // Process monitoring information - extract details for mc and reco neutrinos
    if (m_printToScreen)
        std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    int mcNeutrinoNuance(-1), mcNeutrinoPdg(0), recoNeutrinoPdg(0);
    float mcNeutrinoVtxX(-1.f), mcNeutrinoVtxY(-1.f), mcNeutrinoVtxZ(-1.f);
    float recoNeutrinoVtxX(-1.f), recoNeutrinoVtxY(-1.f), recoNeutrinoVtxZ(-1.f);
    const int nMCNeutrinos(mcNeutrinoList.size()), nRecoNeutrinos(recoNeutrinoList.size()), nMCPrimaries(mcPrimaryList.size());

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        mcNeutrinoPdg = pMCNeutrino->GetParticleId();
        mcNeutrinoVtxX = pMCNeutrino->GetEndpoint().GetX();
        mcNeutrinoVtxY = pMCNeutrino->GetEndpoint().GetY();
        mcNeutrinoVtxZ = pMCNeutrino->GetEndpoint().GetZ();

        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);

        if (pLArMCNeutrino)
            mcNeutrinoNuance = pLArMCNeutrino->GetNuanceCode();

        if (m_printToScreen)
            std::cout << "MCNeutrino, PDG " << mcNeutrinoPdg << ", Nuance " << mcNeutrinoNuance << std::endl;
    }

    for (const ParticleFlowObject *const pPfo : recoNeutrinoList)
    {
        recoNeutrinoPdg = pPfo->GetParticleId();
        const Vertex *const pVertex(pPfo->GetVertexList().empty() ? NULL : *(pPfo->GetVertexList().begin()));
        recoNeutrinoVtxX = pVertex->GetPosition().GetX();
        recoNeutrinoVtxY = pVertex->GetPosition().GetY();
        recoNeutrinoVtxZ = pVertex->GetPosition().GetZ();

        if (m_printToScreen)
            std::cout << "RecoNeutrino, PDG " << recoNeutrinoPdg << std::endl;
    }

    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCNeutrinos", nMCNeutrinos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPdg", mcNeutrinoPdg));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoNuance", mcNeutrinoNuance));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nRecoNeutrinos", nRecoNeutrinos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoPdg", recoNeutrinoPdg));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoVtxX", mcNeutrinoVtxX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoVtxY", mcNeutrinoVtxY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoVtxZ", mcNeutrinoVtxZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxX", recoNeutrinoVtxX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxY", recoNeutrinoVtxY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxZ", recoNeutrinoVtxZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCPrimaries", nMCPrimaries));
    }

    // Process monitoring information - extract details of each mc primary (ordered by number of true hits)
    SimpleMCPrimaryList simpleMCPrimaryList;

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        SimpleMCPrimary simpleMCPrimary;
        simpleMCPrimary.m_pPandoraAddress = pMCPrimary;
        simpleMCPrimary.m_pdgCode = pMCPrimary->GetParticleId();
        simpleMCPrimary.m_energy = pMCPrimary->GetEnergy();
        simpleMCPrimary.m_momentum = pMCPrimary->GetMomentum();
        simpleMCPrimary.m_vertex = pMCPrimary->GetVertex();
        simpleMCPrimary.m_endpoint = pMCPrimary->GetEndpoint();

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

    // Process monitoring information - for each primary, write tree entry containing vector of matched pfos (ordered by number of matched hits)
    int mcPrimaryId(0);
    std::sort(simpleMCPrimaryList.begin(), simpleMCPrimaryList.end(), EventValidationAlgorithm::SortSimpleMCPrimaries);

    for (const SimpleMCPrimary &simpleMCPrimary : simpleMCPrimaryList)
    {
        if (m_printToScreen && (simpleMCPrimary.m_nMCHitsTotal >= m_minHitsToPrintPrimary))
        {
            std::cout << std::endl << "Primary " << mcPrimaryId << ", PDG " << simpleMCPrimary.m_pdgCode << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
                      << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << ")" << std::endl;
        }

        if (m_writeToTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryId", mcPrimaryId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPdg", simpleMCPrimary.m_pdgCode));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsTotal", simpleMCPrimary.m_nMCHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsU", simpleMCPrimary.m_nMCHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsV", simpleMCPrimary.m_nMCHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsW", simpleMCPrimary.m_nMCHitsW));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryE", simpleMCPrimary.m_energy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPX", simpleMCPrimary.m_momentum.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPY", simpleMCPrimary.m_momentum.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPZ", simpleMCPrimary.m_momentum.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxX", simpleMCPrimary.m_vertex.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxY", simpleMCPrimary.m_vertex.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxZ", simpleMCPrimary.m_vertex.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndX", simpleMCPrimary.m_endpoint.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndY", simpleMCPrimary.m_endpoint.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndZ", simpleMCPrimary.m_endpoint.GetZ()));
        }

        // Process monitoring information - first loop over unordered list of matched pfos
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

                // ATTN vertex and end positions/directions currently only filled for track pfos
                const LArTrackPfo *const pLArTrackPfo = dynamic_cast<const LArTrackPfo*>(pMatchedPfo);

                if (pLArTrackPfo)
                {
                    simpleMatchedPfo.m_vertex = pLArTrackPfo->GetVertexPosition();
                    simpleMatchedPfo.m_endpoint = pLArTrackPfo->GetEndPosition();
                    simpleMatchedPfo.m_vertexDirection = pLArTrackPfo->GetVertexDirection();
                    simpleMatchedPfo.m_endDirection = pLArTrackPfo->GetEndDirection();
                }

                simpleMatchedPfoList.push_back(simpleMatchedPfo);
            }
        }

        // Process monitoring information - write the ordered vectors of matched pfo details
        const int mcPrimaryNMatchedPfos(simpleMatchedPfoList.size());
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNMatchedPfos", mcPrimaryNMatchedPfos));

        IntVector pfoIdVector, pfoPdgVector, pfoNHitsTotalVector, pfoNHitsUVector, pfoNHitsVVector, pfoNHitsWVector,
            pfoNMatchedHitsTotalVector, pfoNMatchedHitsUVector, pfoNMatchedHitsVVector, pfoNMatchedHitsWVector;
        FloatVector pfoVtxXVector, pfoVtxYVector, pfoVtxZVector, pfoEndXVector, pfoEndYVector, pfoEndZVector,
            pfoVtxDirXVector, pfoVtxDirYVector, pfoVtxDirZVector, pfoEndDirXVector, pfoEndDirYVector, pfoEndDirZVector;

        std::sort(simpleMatchedPfoList.begin(), simpleMatchedPfoList.end(), EventValidationAlgorithm::SortSimpleMatchedPfos);

        for (const SimpleMatchedPfo simpleMatchedPfo : simpleMatchedPfoList)
        {
            pfoIdVector.push_back(simpleMatchedPfo.m_id);
            pfoPdgVector.push_back(simpleMatchedPfo.m_pdgCode);
            pfoNHitsTotalVector.push_back(simpleMatchedPfo.m_nPfoHitsTotal);
            pfoNHitsUVector.push_back(simpleMatchedPfo.m_nPfoHitsU);
            pfoNHitsVVector.push_back(simpleMatchedPfo.m_nPfoHitsV);
            pfoNHitsWVector.push_back(simpleMatchedPfo.m_nPfoHitsW);
            pfoNMatchedHitsTotalVector.push_back(simpleMatchedPfo.m_nMatchedHitsTotal);
            pfoNMatchedHitsUVector.push_back(simpleMatchedPfo.m_nMatchedHitsU);
            pfoNMatchedHitsVVector.push_back(simpleMatchedPfo.m_nMatchedHitsV);
            pfoNMatchedHitsWVector.push_back(simpleMatchedPfo.m_nMatchedHitsW);

            pfoVtxXVector.push_back(simpleMatchedPfo.m_vertex.GetX());
            pfoVtxYVector.push_back(simpleMatchedPfo.m_vertex.GetY());
            pfoVtxZVector.push_back(simpleMatchedPfo.m_vertex.GetZ());
            pfoEndXVector.push_back(simpleMatchedPfo.m_endpoint.GetX());
            pfoEndYVector.push_back(simpleMatchedPfo.m_endpoint.GetY());
            pfoEndZVector.push_back(simpleMatchedPfo.m_endpoint.GetZ());
            pfoVtxDirXVector.push_back(simpleMatchedPfo.m_vertexDirection.GetX());
            pfoVtxDirYVector.push_back(simpleMatchedPfo.m_vertexDirection.GetY());
            pfoVtxDirZVector.push_back(simpleMatchedPfo.m_vertexDirection.GetZ());
            pfoEndDirXVector.push_back(simpleMatchedPfo.m_endDirection.GetX());
            pfoEndDirYVector.push_back(simpleMatchedPfo.m_endDirection.GetY());
            pfoEndDirZVector.push_back(simpleMatchedPfo.m_endDirection.GetZ());

            if (m_printToScreen && (simpleMCPrimary.m_nMCHitsTotal >= m_minHitsToPrintPrimary))
            {
                std::cout << "-MatchedPfo " << simpleMatchedPfo.m_id << ", PDG " << simpleMatchedPfo.m_pdgCode << ", nMatchedHits " << simpleMatchedPfo.m_nMatchedHitsTotal
                          << " (" << simpleMatchedPfo.m_nMatchedHitsU << ", " << simpleMatchedPfo.m_nMatchedHitsV << ", " << simpleMatchedPfo.m_nMatchedHitsW << ")"
                          << ", nPfoHits " << simpleMatchedPfo.m_nPfoHitsTotal << " (" << simpleMatchedPfo.m_nPfoHitsU << ", " << simpleMatchedPfo.m_nPfoHitsV << ", "
                          << simpleMatchedPfo.m_nPfoHitsW << ")" << std::endl;
            }
        }

        if (m_writeToTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoId", &pfoIdVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoPdg", &pfoPdgVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsTotal", &pfoNHitsTotalVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsU", &pfoNHitsUVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsV", &pfoNHitsVVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNHitsW", &pfoNHitsWVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsTotal", &pfoNMatchedHitsTotalVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsU", &pfoNMatchedHitsUVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsV", &pfoNMatchedHitsVVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoNMatchedHitsW", &pfoNMatchedHitsWVector));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxX", &pfoVtxXVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxY", &pfoVtxYVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxZ", &pfoVtxZVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndX", &pfoEndXVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndY", &pfoEndYVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndZ", &pfoEndZVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxDirX", &pfoVtxDirXVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxDirY", &pfoVtxDirYVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoVtxDirZ", &pfoVtxDirZVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndDirX", &pfoEndDirXVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndDirY", &pfoEndDirYVector));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoEndDirZ", &pfoEndDirZVector));

            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }

        ++mcPrimaryId;
    }

    ++m_eventNumber;
    
    if (m_printToScreen)
        std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
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
    if (lhs.m_nMCHitsTotal != rhs.m_nMCHitsTotal)
        return (lhs.m_nMCHitsTotal > rhs.m_nMCHitsTotal);

    return (lhs.m_energy > rhs.m_energy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::SortSimpleMatchedPfos(const SimpleMatchedPfo &lhs, const SimpleMatchedPfo &rhs)
{
    if (lhs.m_nMatchedHitsTotal != rhs.m_nMatchedHitsTotal)
        return (lhs.m_nMatchedHitsTotal > rhs.m_nMatchedHitsTotal);

    return (lhs.m_nPfoHitsTotal > rhs.m_nPfoHitsTotal);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::SimpleMCPrimary::SimpleMCPrimary() :
    m_pdgCode(0),
    m_nMCHitsTotal(0),  
    m_nMCHitsU(0),
    m_nMCHitsV(0),
    m_nMCHitsW(0),
    m_energy(0.f),
    m_momentum(0.f, 0.f, 0.f),
    m_vertex(-1.f, -1.f, -1.f),
    m_endpoint(-1.f, -1.f, -1.f),
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
    m_vertex(0.f, 0.f, 0.f),
    m_endpoint(0.f, 0.f, 0.f),
    m_vertexDirection(0.f, 0.f, 0.f),
    m_endDirection(0.f, 0.f, 0.f),
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintToScreen", m_printToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsToPrintPrimary", m_minHitsToPrintPrimary));

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
