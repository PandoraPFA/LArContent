/**
 *  @file   LArContent/src/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArPfoHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArMonitoring/EventValidationAlgorithm.h"

#include "LArObjects/LArMCParticle.h"
#include "LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_primaryPfosOnly(true),
    m_collapseToPrimaryPfos(true),
    m_printAllToScreen(false),
    m_printMatchingToScreen(false),
    m_visualizeMatching(false),
    m_writeToTree(true),
    m_matchingMinPrimaryHits(15),
    m_matchingMinSharedHits(5),
    m_vertexVisualizationDeltaR(1.f),
    m_fileIdentifier(0),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::~EventValidationAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::Run()
{
    m_eventNumber++;

    // Input collections
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    
    const VertexList *pTop5VertexList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_top5VertexListName, pTop5VertexList));
    
    const VertexList *pAllVerticesList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_allVerticesListName, pAllVerticesList));

    const PfoList *pPfoList = NULL;
    PfoList inputPfoList((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) ? PfoList(*pPfoList) : PfoList());

    PfoList pfoList;
    LArMonitoringHelper::ExtractTargetPfos(inputPfoList, m_primaryPfosOnly, pfoList);

    // Extract monitoring information
    PfoIdMap pfoIdMap;                                              // pfo -> unique identifier
    this->GetPfoIdMap(pfoList, pfoIdMap);

    MCParticleVector mcNeutrinoList;                                // true neutrinos
    LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoList);

    PfoList recoNeutrinoList;                                       // reco neutrinos
    LArPfoHelper::GetRecoNeutrinos(pPfoList, recoNeutrinoList);

    MCParticleVector mcPrimaryList;                                 // primary mc particles
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryList);

    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;            // [mc particles -> primary mc particle]
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    LArMonitoringHelper::CaloHitToMCMap hitToPrimaryMCMap;          // [hit -> primary mc particle]
    LArMonitoringHelper::MCContributionMap mcToTrueHitListMap;      // [primary mc particle -> true hit list]
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

    LArMonitoringHelper::CaloHitToPfoMap hitToPfoMap;               // [hit -> pfo]
    LArMonitoringHelper::PfoContributionMap pfoToHitListMap;        // [pfo -> reco hit list]
    LArMonitoringHelper::GetPfoToCaloHitMatches(pCaloHitList, pfoList, m_collapseToPrimaryPfos, hitToPfoMap, pfoToHitListMap);

    LArMonitoringHelper::MCToPfoMap mcToBestPfoMap;                 // [mc particle -> best matched pfo]
    LArMonitoringHelper::MCContributionMap mcToBestPfoHitsMap;      // [mc particle -> list of hits included in best pfo]
    LArMonitoringHelper::MCToPfoMatchingMap mcToFullPfoMatchingMap; // [mc particle -> all matched pfos (and matched hits)]
    LArMonitoringHelper::GetMCParticleToPfoMatches(pCaloHitList, pfoToHitListMap, hitToPrimaryMCMap, mcToBestPfoMap, mcToBestPfoHitsMap, mcToFullPfoMatchingMap);

    SimpleMCPrimaryList simpleMCPrimaryList;
    this->GetSimpleMCPrimaryList(mcPrimaryList, mcToTrueHitListMap, mcToFullPfoMatchingMap, simpleMCPrimaryList);

    MCPrimaryMatchingMap mcPrimaryMatchingMap;
    this->GetMCPrimaryMatchingMap(simpleMCPrimaryList, pfoIdMap, mcToFullPfoMatchingMap, pfoToHitListMap, mcPrimaryMatchingMap);

    if (m_printAllToScreen)
        this->PrintAllOutput(mcNeutrinoList, recoNeutrinoList, mcPrimaryMatchingMap);

    if (m_writeToTree)
        this->WriteAllOutput(mcNeutrinoList, recoNeutrinoList, mcPrimaryMatchingMap, pTop5VertexList, pAllVerticesList);

    if (m_printMatchingToScreen || m_visualizeMatching)
    {
        MatchingDetailsMap matchingDetailsMap;
        this->PerformMatching(mcPrimaryMatchingMap, matchingDetailsMap);

        if (m_printMatchingToScreen)
            this->PrintMatchingOutput(mcPrimaryMatchingMap, matchingDetailsMap);

        if (m_visualizeMatching)
            this->VisualizeMatchingOutput(mcNeutrinoList, recoNeutrinoList, mcPrimaryMatchingMap, matchingDetailsMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetSimpleMCPrimaryList(const MCParticleVector &mcPrimaryList, const LArMonitoringHelper::MCContributionMap &mcToTrueHitListMap,
    const LArMonitoringHelper::MCToPfoMatchingMap &mcToFullPfoMatchingMap, SimpleMCPrimaryList &simpleMCPrimaryList) const
{
    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        SimpleMCPrimary simpleMCPrimary;
        // ATTN simpleMCPrimary.m_id assigned later, after sorting
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

    std::sort(simpleMCPrimaryList.begin(), simpleMCPrimaryList.end(), EventValidationAlgorithm::SortSimpleMCPrimaries);

    int mcPrimaryId(0);
    for (SimpleMCPrimary &simpleMCPrimary : simpleMCPrimaryList)
        simpleMCPrimary.m_id = mcPrimaryId++;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetMCPrimaryMatchingMap(const SimpleMCPrimaryList &simpleMCPrimaryList, const PfoIdMap &pfoIdMap,
    const LArMonitoringHelper::MCToPfoMatchingMap &mcToFullPfoMatchingMap, const LArMonitoringHelper::PfoContributionMap &pfoToHitListMap,
    MCPrimaryMatchingMap &mcPrimaryMatchingMap) const
{
    for (const SimpleMCPrimary &simpleMCPrimary : simpleMCPrimaryList)
    {
        // First loop over unordered list of matched pfos
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

                // ATTN Assume pfos have either zero or one parents. Ignore parent neutrino.
                const ParticleFlowObject *const pParentPfo(pMatchedPfo->GetParentPfoList().empty() ? NULL : *(pMatchedPfo->GetParentPfoList().begin()));

                if (pParentPfo && !LArPfoHelper::IsNeutrino(pParentPfo))
                    simpleMatchedPfo.m_parentId = pfoIdMap.at(pParentPfo);

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

        // Store the ordered vectors of matched pfo details
        std::sort(simpleMatchedPfoList.begin(), simpleMatchedPfoList.end(), EventValidationAlgorithm::SortSimpleMatchedPfos);

        if (!mcPrimaryMatchingMap.insert(MCPrimaryMatchingMap::value_type(simpleMCPrimary, simpleMatchedPfoList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintAllOutput(const MCParticleVector &mcNeutrinoList, const PfoList &recoNeutrinoList,
    const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const
{
    std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);
        std::cout << "MCNeutrino, PDG " << pMCNeutrino->GetParticleId() << ", Nuance " << (pLArMCNeutrino ? pLArMCNeutrino->GetNuanceCode() : -1) << std::endl;
    }

    for (const ParticleFlowObject *const pPfo : recoNeutrinoList)
    {
        std::cout << "RecoNeutrino, PDG " << pPfo->GetParticleId() << std::endl;

        if ((1 == pPfo->GetVertexList().size()) && (1 == mcNeutrinoList.size()))
            std::cout << "VtxOffset" << ((*(pPfo->GetVertexList().begin()))->GetPosition() - (*(mcNeutrinoList.begin()))->GetEndpoint()) << std::endl;
    }

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        std::cout << std::endl << "Primary " << simpleMCPrimary.m_id << ", PDG " << simpleMCPrimary.m_pdgCode << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
            << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << ")" << std::endl;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            std::cout << "-MatchedPfo " << simpleMatchedPfo.m_id;

            if (simpleMatchedPfo.m_parentId >= 0)
                std::cout << ", ParentPfo " << simpleMatchedPfo.m_parentId;

            std::cout << ", PDG " << simpleMatchedPfo.m_pdgCode << ", nMatchedHits " << simpleMatchedPfo.m_nMatchedHitsTotal
                << " (" << simpleMatchedPfo.m_nMatchedHitsU << ", " << simpleMatchedPfo.m_nMatchedHitsV << ", " << simpleMatchedPfo.m_nMatchedHitsW << ")"
                << ", nPfoHits " << simpleMatchedPfo.m_nPfoHitsTotal << " (" << simpleMatchedPfo.m_nPfoHitsU << ", " << simpleMatchedPfo.m_nPfoHitsV << ", "
                << simpleMatchedPfo.m_nPfoHitsW << ")" << std::endl;
        }
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::WriteAllOutput(const MCParticleVector &mcNeutrinoList, const PfoList &recoNeutrinoList,
    const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const VertexList* pTop5VertexList, const VertexList* pAllVerticesList) const
{
#ifdef MONITORING
    int mcNeutrinoNuance(-1), mcNeutrinoPdg(0), recoNeutrinoPdg(0);
    float mcNeutrinoVtxX(-1.f), mcNeutrinoVtxY(-1.f), mcNeutrinoVtxZ(-1.f);
    float recoNeutrinoVtxX(-1.f), recoNeutrinoVtxY(-1.f), recoNeutrinoVtxZ(-1.f);
    float mcNeutrinoE(0.f), mcNeutrinoPX(0.f), mcNeutrinoPY(0.f), mcNeutrinoPZ(0.f);
    const int nMCNeutrinos(mcNeutrinoList.size()), nRecoNeutrinos(recoNeutrinoList.size()), nMCPrimaries(mcPrimaryMatchingMap.size());

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        mcNeutrinoPdg = pMCNeutrino->GetParticleId();
        mcNeutrinoVtxX = pMCNeutrino->GetEndpoint().GetX();
        mcNeutrinoVtxY = pMCNeutrino->GetEndpoint().GetY();
        mcNeutrinoVtxZ = pMCNeutrino->GetEndpoint().GetZ();

        mcNeutrinoE = pMCNeutrino->GetEnergy();
        mcNeutrinoPX = pMCNeutrino->GetMomentum().GetX();
        mcNeutrinoPY = pMCNeutrino->GetMomentum().GetY();
        mcNeutrinoPZ = pMCNeutrino->GetMomentum().GetZ();

        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);

        if (pLArMCNeutrino)
            mcNeutrinoNuance = pLArMCNeutrino->GetNuanceCode();
    }

    for (const ParticleFlowObject *const pPfo : recoNeutrinoList)
    {
        recoNeutrinoPdg = pPfo->GetParticleId();
        const Vertex *const pVertex(pPfo->GetVertexList().empty() ? NULL : *(pPfo->GetVertexList().begin()));
        recoNeutrinoVtxX = pVertex->GetPosition().GetX();
        recoNeutrinoVtxY = pVertex->GetPosition().GetY();
        recoNeutrinoVtxZ = pVertex->GetPosition().GetZ();
    }
    
    //---------------------------------------------------------TOP 5--------------------------------------------------------------------
    if (recoNeutrinoList.size() == 1 && mcNeutrinoList.size() == 1)
    {
        float top5VertexOffset(-1.f), bestVertexOffset(-1.f), top5VertexOffsetX(-1.f), top5VertexOffsetY(-1.f), top5VertexOffsetZ(-1.f), bestVertexOffsetX(-1.f),  bestVertexOffsetY(-1.f),  bestVertexOffsetZ(-1.f);
        
        std::vector<float> top5VerticesDR;
        std::vector<float> allVerticesDR;
        
        CartesianVector mcNeutrinoVertexPosition((*mcNeutrinoList.front()).GetVertex());
        
        //const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), mcNeutrinoVertexPosition, TPC_VIEW_U));
        //const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), mcNeutrinoVertexPosition, TPC_VIEW_V));
        //const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), mcNeutrinoVertexPosition, TPC_VIEW_W));
        //
        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionU, "Target Vertex", RED, 1));
        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionV, "Target Vertex", RED, 1));
        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionW, "Target Vertex", RED, 1));
        //
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        
        if (!pTop5VertexList->empty())
        {
            for (const Vertex *const pVertex : (*pTop5VertexList))
            {
                float vertexDR((pVertex->GetPosition() - mcNeutrinoVertexPosition).GetMagnitude());
                top5VerticesDR.push_back(vertexDR);
            }
            
            top5VertexOffset = (*std::min_element(top5VerticesDR.begin(), top5VerticesDR.end()));
            
            for (const Vertex *const pVertex : (*pTop5VertexList))
            {
                float vertexDR((pVertex->GetPosition() - mcNeutrinoVertexPosition).GetMagnitude());
                if (vertexDR == top5VertexOffset)
                {
                    top5VertexOffsetX = (mcNeutrinoVertexPosition.GetX() - pVertex->GetPosition().GetX());
                    top5VertexOffsetY = (mcNeutrinoVertexPosition.GetY() - pVertex->GetPosition().GetY());
                    top5VertexOffsetZ = (mcNeutrinoVertexPosition.GetZ() - pVertex->GetPosition().GetZ());
                }
            }
            
            std::cout << "Top 5 vertex DR: " << top5VertexOffset << " with DX: " << top5VertexOffsetX << " DY: " << top5VertexOffsetY << " DZ: " << top5VertexOffsetZ << std::endl;
        }
        
        if (!pAllVerticesList->empty())
        {
            std::cout << "All vertex list contains " << pAllVerticesList->size() << " entries." << std::endl;
            
            for (const Vertex *const pVertex : (*pAllVerticesList))
            {
                float vertexDR((pVertex->GetPosition() - mcNeutrinoVertexPosition).GetMagnitude());
                allVerticesDR.push_back(vertexDR);
            }
            
            bestVertexOffset = (*std::min_element(allVerticesDR.begin(), allVerticesDR.end()));
            
            for (const Vertex *const pVertex : (*pAllVerticesList))
            {
                float vertexDR((pVertex->GetPosition() - mcNeutrinoVertexPosition).GetMagnitude());
                if (vertexDR == bestVertexOffset)
                {
                    bestVertexOffsetX = (mcNeutrinoVertexPosition.GetX() - pVertex->GetPosition().GetX());
                    bestVertexOffsetY = (mcNeutrinoVertexPosition.GetY() - pVertex->GetPosition().GetY());
                    bestVertexOffsetZ = (mcNeutrinoVertexPosition.GetZ() - pVertex->GetPosition().GetZ());
                }
            }
            
            std::cout << "Best possible vertex DR: " << bestVertexOffset << " with DX: " << bestVertexOffsetX << " DY: " << bestVertexOffsetY << " DZ: " << bestVertexOffsetZ << std::endl;
        }
        
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "top5VertexOffset", top5VertexOffset));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestVertexOffset", bestVertexOffset));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "top5VertexOffsetX", top5VertexOffsetX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "top5VertexOffsetY", top5VertexOffsetY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "top5VertexOffsetZ", top5VertexOffsetZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestVertexOffsetX", bestVertexOffsetX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestVertexOffsetY", bestVertexOffsetY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestVertexOffsetZ", bestVertexOffsetZ));
    }
//---------------------------------------------------------TOP 5--------------------------------------------------------------------

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
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoE", mcNeutrinoE));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPX", mcNeutrinoPX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPY", mcNeutrinoPY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoPZ", mcNeutrinoPZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxX", recoNeutrinoVtxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxY", recoNeutrinoVtxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNeutrinoVtxZ", recoNeutrinoVtxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCPrimaries", nMCPrimaries));

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryId", simpleMCPrimary.m_id));
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

        const int mcPrimaryNMatchedPfos(mapValue.second.size());
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNMatchedPfos", mcPrimaryNMatchedPfos));

        IntVector pfoIdVector, pfoParentIdVector, pfoPdgVector, pfoNHitsTotalVector, pfoNHitsUVector, pfoNHitsVVector, pfoNHitsWVector,
            pfoNMatchedHitsTotalVector, pfoNMatchedHitsUVector, pfoNMatchedHitsVVector, pfoNMatchedHitsWVector;
        FloatVector pfoVtxXVector, pfoVtxYVector, pfoVtxZVector, pfoEndXVector, pfoEndYVector, pfoEndZVector,
            pfoVtxDirXVector, pfoVtxDirYVector, pfoVtxDirZVector, pfoEndDirXVector, pfoEndDirYVector, pfoEndDirZVector;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            pfoIdVector.push_back(simpleMatchedPfo.m_id);
            pfoParentIdVector.push_back(simpleMatchedPfo.m_parentId);
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
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoId", &pfoIdVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "matchedPfoParentId", &pfoParentIdVector));
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
#else
    std::cout << "Monitoring functionality unavailable. nMCNeutrinos " << mcNeutrinoList.size()
              << ", nRecoNeutrinos" << recoNeutrinoList.size() << ", nMCPrimaries " << mcPrimaryMatchingMap.size() << std::endl;
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PerformMatching(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, MatchingDetailsMap &matchingDetailsMap) const
{
    // Get best matches, one-by-one, until no more strong matches possible
    IntSet usedMCIds, usedPfoIds;
    while (GetStrongestPfoMatch(mcPrimaryMatchingMap, usedMCIds, usedPfoIds, matchingDetailsMap)) {}

    // Assign any remaining pfos to primaries, based on number of matched hits
    GetRemainingPfoMatches(mcPrimaryMatchingMap, usedPfoIds, matchingDetailsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::GetStrongestPfoMatch(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, IntSet &usedMCIds, IntSet &usedPfoIds,
    MatchingDetailsMap &matchingDetailsMap) const
{
    int bestPfoMatchId(-1);
    MatchingDetails bestMatchingDetails;

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        if (simpleMCPrimary.m_nMCHitsTotal < m_matchingMinPrimaryHits)
            continue;

        if (usedMCIds.count(simpleMCPrimary.m_id))
            continue;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (usedPfoIds.count(simpleMatchedPfo.m_id) || (simpleMatchedPfo.m_nMatchedHitsTotal < m_matchingMinSharedHits))
                continue;

            if (simpleMatchedPfo.m_nMatchedHitsTotal > bestMatchingDetails.m_nMatchedHits)
            {
                bestPfoMatchId = simpleMatchedPfo.m_id;
                bestMatchingDetails.m_matchedPrimaryId = simpleMCPrimary.m_id;
                bestMatchingDetails.m_nMatchedHits = simpleMatchedPfo.m_nMatchedHitsTotal;
                bestMatchingDetails.m_completeness = static_cast<float>(simpleMatchedPfo.m_nMatchedHitsTotal) / static_cast<float>(simpleMCPrimary.m_nMCHitsTotal);
            }
        }
    }

    if (bestPfoMatchId > -1)
    {
        matchingDetailsMap[bestPfoMatchId] = bestMatchingDetails;
        usedMCIds.insert(bestMatchingDetails.m_matchedPrimaryId);
        usedPfoIds.insert(bestPfoMatchId);
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetRemainingPfoMatches(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const IntSet &usedPfoIds,
    MatchingDetailsMap &matchingDetailsMap) const
{
    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        if (simpleMCPrimary.m_nMCHitsTotal < m_matchingMinPrimaryHits)
            continue;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (usedPfoIds.count(simpleMatchedPfo.m_id) || (simpleMatchedPfo.m_nMatchedHitsTotal < m_matchingMinSharedHits))
                continue;

            MatchingDetails &matchingDetails(matchingDetailsMap[simpleMatchedPfo.m_id]);

            if (simpleMatchedPfo.m_nMatchedHitsTotal > matchingDetails.m_nMatchedHits)
            {
                matchingDetails.m_matchedPrimaryId = simpleMCPrimary.m_id;
                matchingDetails.m_nMatchedHits = simpleMatchedPfo.m_nMatchedHitsTotal;
                matchingDetails.m_completeness = static_cast<float>(simpleMatchedPfo.m_nMatchedHitsTotal) / static_cast<float>(simpleMCPrimary.m_nMCHitsTotal);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintMatchingOutput(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const
{
    std::cout << "---PROCESSED-MATCHING-OUTPUT--------------------------------------------------------------------" << std::endl;
    bool isCorrect(true);

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        if (simpleMCPrimary.m_nMCHitsTotal < m_matchingMinPrimaryHits)
            continue;

        std::cout << std::endl << "Primary " << simpleMCPrimary.m_id << ", PDG " << simpleMCPrimary.m_pdgCode << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
            << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << ")" << std::endl;

        unsigned int nMatches(0);

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
            {
                std::cout << "-MatchedPfo " << simpleMatchedPfo.m_id;
                ++nMatches;

                if (simpleMatchedPfo.m_parentId >= 0)
                    std::cout << ", ParentPfo " << simpleMatchedPfo.m_parentId;

                std::cout << ", PDG " << simpleMatchedPfo.m_pdgCode << ", nMatchedHits " << simpleMatchedPfo.m_nMatchedHitsTotal
                    << " (" << simpleMatchedPfo.m_nMatchedHitsU << ", " << simpleMatchedPfo.m_nMatchedHitsV << ", " << simpleMatchedPfo.m_nMatchedHitsW << ")"
                    << ", nPfoHits " << simpleMatchedPfo.m_nPfoHitsTotal << " (" << simpleMatchedPfo.m_nPfoHitsU << ", " << simpleMatchedPfo.m_nPfoHitsV << ", "
                    << simpleMatchedPfo.m_nPfoHitsW << ")" << std::endl;
            }
        }

        if (1 != nMatches)
            isCorrect = false;
    }

    std::cout << std::endl << "Is correct? " << isCorrect << std::endl;
    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::VisualizeMatchingOutput(const MCParticleVector &mcNeutrinoList, const PfoList &recoNeutrinoList,
    const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const
{
    this->VisualizeVertexMatches(mcNeutrinoList, recoNeutrinoList);
    unsigned int displayIndex(0);

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        if (simpleMCPrimary.m_nMCHitsTotal < m_matchingMinPrimaryHits)
            continue;

        const std::string displayString(TypeToString(displayIndex));
        PfoList primaryMatchedPfos, matchedPfos;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
            {
                primaryMatchedPfos.insert(simpleMatchedPfo.m_pPandoraAddress);
                LArPfoHelper::GetAllDownstreamPfos(simpleMatchedPfo.m_pPandoraAddress, matchedPfos);
            }
        }

#ifdef MONITORING
        const std::string prefix((1 == primaryMatchedPfos.size()) ? "Matched" : "Split");
        const Color color((1 == primaryMatchedPfos.size()) ? GREEN : RED);

        for (const ParticleFlowObject *const pPfo : matchedPfos)
        {
            ClusterList clusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
            //PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &matchedPfos, prefix + "Pfo_" + displayString, color, true, false));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, prefix + "PfoClusters_" + displayString, color));
            PANDORA_MONITORING_API(VisualizeVertices(this->GetPandora(), &(pPfo->GetVertexList()), prefix + "PfoVertices_" + displayString, color));
        }

        if (matchedPfos.empty())
        {
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &simpleMCPrimary.m_vertex, "MissingPrimaryVtx_" + displayString, RED, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &simpleMCPrimary.m_endpoint, "MissingPrimaryEnd_" + displayString, RED, 1));
        }
#endif
        ++displayIndex;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::VisualizeVertexMatches(const MCParticleVector &mcNeutrinoList, const PfoList &recoNeutrinoList) const
{
#ifdef MONITORING
    const Vertex *pBestVertex(nullptr);
    float closestDistance(m_vertexVisualizationDeltaR);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        const CartesianVector &mcVertexPosition(pMCNeutrino->GetEndpoint());

        for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinoList)
        {
            const Vertex *const pVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
            const float distance((mcVertexPosition - pVertex->GetPosition()).GetMagnitude());

            if (distance < closestDistance)
            {
                closestDistance = distance;
                pBestVertex = pVertex;
            }
        }
    }

    const bool isGood((1 == mcNeutrinoList.size()) && (1 == recoNeutrinoList.size()) && pBestVertex);
    const Color color(isGood ? GRAY : BLUE);

    unsigned int displayIndex(0);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoList)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &(pMCNeutrino->GetEndpoint()), "MCNeutrinoVertex_" + TypeToString(displayIndex++), color, 1));
    }

    displayIndex = 0;

    for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinoList)
    {
        const Vertex *const pThisVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));

        const bool isThisGood(isGood || ((1 == mcNeutrinoList.size()) && (pThisVertex == pBestVertex)));
        const std::string thisPrefix(isThisGood ? "Good" : "Displaced");
        const Color thisColor(isThisGood ? CYAN : VIOLET);

        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &(pThisVertex->GetPosition()), thisPrefix + "RecoNeutrinoVertex_" + TypeToString(displayIndex++), thisColor, 1));
    }
#else
    std::cout << "Monitoring functionality unavailable. nMCNeutrinos " << mcNeutrinoList.size() << ", nRecoNeutrinos" << recoNeutrinoList.size() << std::endl;
#endif
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

    if (lhs.m_nPfoHitsTotal != rhs.m_nPfoHitsTotal)
        return (lhs.m_nPfoHitsTotal > rhs.m_nPfoHitsTotal);

    return (lhs.m_id < rhs.m_id);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::SimpleMCPrimary::SimpleMCPrimary() :
    m_id(-1),
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

bool EventValidationAlgorithm::SimpleMCPrimary::operator<(const SimpleMCPrimary &rhs) const
{
    if (this == &rhs)
        return false;

    return (m_id < rhs.m_id);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::SimpleMatchedPfo::SimpleMatchedPfo() :
    m_id(-1),
    m_parentId(-1),
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

EventValidationAlgorithm::MatchingDetails::MatchingDetails() :
    m_matchedPrimaryId(-1),
    m_nMatchedHits(0),
    m_completeness(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Top5VertexListName", m_top5VertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "AllVerticesListName", m_allVerticesListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrimaryPfosOnly", m_primaryPfosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CollapseToPrimaryPfos", m_collapseToPrimaryPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeMatching", m_visualizeMatching));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPrimaryHits", m_matchingMinPrimaryHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinSharedHits", m_matchingMinSharedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexVisualizationDeltaR", m_vertexVisualizationDeltaR));

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
