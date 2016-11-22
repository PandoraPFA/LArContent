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
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_integrateOverSlices(false),
    m_neutrinoInducedOnly(false),
    m_primaryPfosOnly(true),
    m_collapseToPrimaryPfos(true),
    m_minHitSharingFraction(0.9f),
    m_maxPhotonPropagation(2.5f),
    m_printAllToScreen(false),
    m_printMatchingToScreen(false),
    m_visualizeMatching(false),
    m_visualizeVertices(false),
    m_visualizeRemnants(false),
    m_visualizeGaps(false),
    m_writeToTree(true),
    m_matchingMinPrimaryHits(15),
    m_matchingMinHitsForGoodView(5),
    m_matchingMinPrimaryGoodViews(2),
    m_useSmallPrimaries(true),
    m_matchingMinSharedHits(5),
    m_matchingMinCompleteness(0.1f),
    m_matchingMinPurity(0.5f),
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
    ++m_eventNumber;

    // Extract input collections
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    // Obtain vector: reco neutrino(s)
    PfoList allRecoNeutrinoList;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, allRecoNeutrinoList);
    allRecoNeutrinoList.sort(EventValidationAlgorithm::SortRecoNeutrinos);

    // Unless integrating over all slices, take only first reco neutrino after sorting. May need a more careful selection here cosmics bkg.
    const PfoList recoNeutrinoList(m_integrateOverSlices ? allRecoNeutrinoList : !allRecoNeutrinoList.empty() ? PfoList(1, allRecoNeutrinoList.front()) : PfoList());
    const PfoVector recoNeutrinoVector(recoNeutrinoList.begin(), recoNeutrinoList.end());

    // Obtain vector: target pfos
    PfoList pfoList; // TODO config
    LArMonitoringHelper::ExtractTargetPfos(recoNeutrinoList, m_primaryPfosOnly, pfoList);

    // Obtain map: pfo -> unique identifier
    PfoIdMap pfoIdMap;
    this->GetPfoIdMap(pfoList, pfoIdMap);

    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoVector);

    // Obtain vector: primary mc particles
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);

    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList selectedCaloHitList;
    this->SelectCaloHits(pCaloHitList, mcToPrimaryMCMap, selectedCaloHitList);

    // Obtain maps: [hit -> primary mc particle], [primary mc particle -> list of hits]
    LArMonitoringHelper::CaloHitToMCMap hitToPrimaryMCMap;
    LArMonitoringHelper::MCContributionMap mcToTrueHitListMap;
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

    // Obtain maps: [hit -> pfo], [pfo -> list of hits]
    LArMonitoringHelper::CaloHitToPfoMap hitToPfoMap;
    LArMonitoringHelper::PfoContributionMap pfoToHitListMap;
    LArMonitoringHelper::GetPfoToCaloHitMatches(&selectedCaloHitList, pfoList, m_collapseToPrimaryPfos, hitToPfoMap, pfoToHitListMap);

    // Obtain maps: [mc particle -> best matched pfo], [mc particle -> list of hits included in best pfo], [mc particle -> all matched pfos (and matched hits)]
    LArMonitoringHelper::MCToPfoMap mcToBestPfoMap;
    LArMonitoringHelper::MCContributionMap mcToBestPfoHitsMap;
    LArMonitoringHelper::MCToPfoMatchingMap mcToFullPfoMatchingMap;
    LArMonitoringHelper::GetMCParticleToPfoMatches(&selectedCaloHitList, pfoToHitListMap, hitToPrimaryMCMap, mcToBestPfoMap, mcToBestPfoHitsMap, mcToFullPfoMatchingMap);

    // Remove shared hits where target particle deposits below threshold energy fraction
    CaloHitList goodCaloHitList;
    this->SelectGoodCaloHits(&selectedCaloHitList, mcToPrimaryMCMap, goodCaloHitList);

    // Obtain maps: [good hit -> primary mc particle], [primary mc particle -> list of good hits]
    LArMonitoringHelper::CaloHitToMCMap goodHitToPrimaryMCMap;
    LArMonitoringHelper::MCContributionMap mcToGoodTrueHitListMap;
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(&goodCaloHitList, mcToPrimaryMCMap, goodHitToPrimaryMCMap, mcToGoodTrueHitListMap);

    // Obtain vector: simple mc primaries
    SimpleMCPrimaryList simpleMCPrimaryList;
    this->GetSimpleMCPrimaryList(mcPrimaryVector, mcToTrueHitListMap, mcToGoodTrueHitListMap, mcToFullPfoMatchingMap, simpleMCPrimaryList);

    // Obtain map: [simple mc primary -> list of simple matched pfos]
    MCPrimaryMatchingMap mcPrimaryMatchingMap;
    this->GetMCPrimaryMatchingMap(simpleMCPrimaryList, pfoIdMap, mcToFullPfoMatchingMap, pfoToHitListMap, mcPrimaryMatchingMap);

    // Print raw matching information to terminal
    if (m_printAllToScreen)
        this->PrintAllOutput(mcNeutrinoVector, recoNeutrinoVector, mcPrimaryMatchingMap);

    // Write raw matching information to root file
    if (m_writeToTree)
        this->WriteAllOutput(mcNeutrinoVector, recoNeutrinoVector, mcPrimaryMatchingMap);

    if (m_printMatchingToScreen || m_visualizeMatching)
    {
        // Obtain map: [simple mc primary -> interpreted list of simple matched pfos]
        MatchingDetailsMap matchingDetailsMap;
        this->PerformMatching(mcPrimaryMatchingMap, matchingDetailsMap);

        // Print interpreted matching information to terminal
        if (m_printMatchingToScreen)
            this->PrintMatchingOutput(mcPrimaryMatchingMap, matchingDetailsMap);
#ifdef MONITORING
        // Visualize interpreted matching information
        if (m_visualizeMatching)
            this->VisualizeMatchingOutput(mcNeutrinoVector, recoNeutrinoVector, mcPrimaryMatchingMap, matchingDetailsMap);
#endif
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedCaloHitList) const
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

            if (mcToPrimaryMCMap.end() == mcIter)
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;

            if (this->PassMCParticleChecks(pPrimaryParticle, pPrimaryParticle, pHitParticle))
                selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedGoodCaloHitList) const
{
    for (const CaloHit *const pCaloHit : *pSelectedCaloHitList)
    {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        MCParticleWeightMap primaryWeightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pMCParticle);

            if (mcToPrimaryMCMap.end() != mcIter)
                primaryWeightMap[mcIter->second] += weight;
        }

        MCParticleVector mcPrimaryVector;
        for (const auto &mapEntry : primaryWeightMap) mcPrimaryVector.push_back(mapEntry.first);
        std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), PointerLessThan<MCParticle>());

        const MCParticle *pBestPrimaryParticle(nullptr);
        float bestPrimaryWeight(0.f), primaryWeightSum(0.f);

        for (const MCParticle *const pPrimaryMCParticle : mcPrimaryVector)
        {
            const float primaryWeight(primaryWeightMap.at(pPrimaryMCParticle));
            primaryWeightSum += primaryWeight;

            if (primaryWeight > bestPrimaryWeight)
            {
                bestPrimaryWeight = primaryWeight;
                pBestPrimaryParticle = pPrimaryMCParticle;
            }
        }

        if (!pBestPrimaryParticle || (primaryWeightSum < std::numeric_limits<float>::epsilon()) || ((bestPrimaryWeight / primaryWeightSum) < m_minHitSharingFraction))
            continue;

        selectedGoodCaloHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
    const MCParticle *const pHitMCParticle) const
{
    if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
        return false;

    if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != pOriginalPrimary->GetParticleId()) && (E_MINUS != std::abs(pOriginalPrimary->GetParticleId())))
    {
        if ((pThisMCParticle->GetEndpoint() - pThisMCParticle->GetVertex()).GetMagnitude() > m_maxPhotonPropagation)
            return false;
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (this->PassMCParticleChecks(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetSimpleMCPrimaryList(const MCParticleVector &mcPrimaryVector, const LArMonitoringHelper::MCContributionMap &mcToTrueHitListMap,
    const LArMonitoringHelper::MCContributionMap &mcToGoodTrueHitListMap, const LArMonitoringHelper::MCToPfoMatchingMap &mcToFullPfoMatchingMap,
    SimpleMCPrimaryList &simpleMCPrimaryList) const
{
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (m_neutrinoInducedOnly && !LArMCParticleHelper::IsNeutrinoInduced(pMCPrimary))
            continue;

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

        LArMonitoringHelper::MCContributionMap::const_iterator goodTrueHitsIter = mcToGoodTrueHitListMap.find(pMCPrimary);

        if (mcToGoodTrueHitListMap.end() != goodTrueHitsIter)
        {
            const CaloHitList &caloHitList(goodTrueHitsIter->second);
            simpleMCPrimary.m_nGoodMCHitsTotal = caloHitList.size();
            simpleMCPrimary.m_nGoodMCHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList);
            simpleMCPrimary.m_nGoodMCHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList);
            simpleMCPrimary.m_nGoodMCHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList);
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
                const ParticleFlowObject *const pParentPfo(pMatchedPfo->GetParentPfoList().empty() ? nullptr : *(pMatchedPfo->GetParentPfoList().begin()));

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

void EventValidationAlgorithm::PrintAllOutput(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector,
    const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const
{
    std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);
        std::cout << "MCNeutrino, PDG " << pMCNeutrino->GetParticleId() << ", Nuance " << (pLArMCNeutrino ? pLArMCNeutrino->GetNuanceCode() : -1) << std::endl;
    }

    for (const ParticleFlowObject *const pPfo : recoNeutrinoVector)
    {
        std::cout << "RecoNeutrino, PDG " << pPfo->GetParticleId() << std::endl;

        if ((1 == pPfo->GetVertexList().size()) && (1 == mcNeutrinoVector.size()))
            std::cout << "VtxOffset" << ((*(pPfo->GetVertexList().begin()))->GetPosition() - mcNeutrinoVector.front()->GetEndpoint()) << std::endl;
    }

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);

        std::cout << std::endl << "Primary " << simpleMCPrimary.m_id << ", PDG " << simpleMCPrimary.m_pdgCode << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
            << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << "),"
            << " [nGood " << simpleMCPrimary.m_nGoodMCHitsTotal << " (" << simpleMCPrimary.m_nGoodMCHitsU << ", " << simpleMCPrimary.m_nGoodMCHitsV
            << ", " << simpleMCPrimary.m_nGoodMCHitsW << ")]" << std::endl;

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

void EventValidationAlgorithm::WriteAllOutput(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector,
    const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const
{
#ifdef MONITORING
    int mcNeutrinoNuance(-1), mcNeutrinoPdg(0), recoNeutrinoPdg(0);
    float mcNeutrinoVtxX(-1.f), mcNeutrinoVtxY(-1.f), mcNeutrinoVtxZ(-1.f);
    float recoNeutrinoVtxX(-1.f), recoNeutrinoVtxY(-1.f), recoNeutrinoVtxZ(-1.f);
    float mcNeutrinoE(0.f), mcNeutrinoPX(0.f), mcNeutrinoPY(0.f), mcNeutrinoPZ(0.f);
    const int nMCNeutrinos(mcNeutrinoVector.size()), nRecoNeutrinos(recoNeutrinoVector.size()), nMCPrimaries(mcPrimaryMatchingMap.size());

    if (!mcNeutrinoVector.empty())
    {
        const MCParticle *const pMCNeutrino = mcNeutrinoVector.front();

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

    if (!recoNeutrinoVector.empty())
    {
        const ParticleFlowObject *const pPfo = recoNeutrinoVector.front();

        recoNeutrinoPdg = pPfo->GetParticleId();
        const Vertex *const pVertex(pPfo->GetVertexList().empty() ? nullptr : *(pPfo->GetVertexList().begin()));

        if (pVertex)
        {
            recoNeutrinoVtxX = pVertex->GetPosition().GetX();
            recoNeutrinoVtxY = pVertex->GetPosition().GetY();
            recoNeutrinoVtxZ = pVertex->GetPosition().GetZ();
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
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
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsTotal", simpleMCPrimary.m_nGoodMCHitsTotal));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsU", simpleMCPrimary.m_nGoodMCHitsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsV", simpleMCPrimary.m_nGoodMCHitsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNGoodHitsW", simpleMCPrimary.m_nGoodMCHitsW));

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
    std::cout << "Monitoring functionality unavailable. nMCNeutrinos " << mcNeutrinoVector.size()
              << ", nRecoNeutrinos" << recoNeutrinoVector.size() << ", nMCPrimaries " << mcPrimaryMatchingMap.size() << std::endl;
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

        if (!m_useSmallPrimaries && !this->IsGoodMCPrimary(simpleMCPrimary))
            continue;

        if (usedMCIds.count(simpleMCPrimary.m_id))
            continue;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (usedPfoIds.count(simpleMatchedPfo.m_id))
                continue;

            if (!this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo))
                continue;

            if (simpleMatchedPfo.m_nMatchedHitsTotal > bestMatchingDetails.m_nMatchedHits)
            {
                bestPfoMatchId = simpleMatchedPfo.m_id;
                bestMatchingDetails.m_matchedPrimaryId = simpleMCPrimary.m_id;
                bestMatchingDetails.m_nMatchedHits = simpleMatchedPfo.m_nMatchedHitsTotal;
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

        if (!m_useSmallPrimaries && !this->IsGoodMCPrimary(simpleMCPrimary))
            continue;

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (usedPfoIds.count(simpleMatchedPfo.m_id))
                continue;

            MatchingDetails &matchingDetails(matchingDetailsMap[simpleMatchedPfo.m_id]);

            if (simpleMatchedPfo.m_nMatchedHitsTotal > matchingDetails.m_nMatchedHits)
            {
                matchingDetails.m_matchedPrimaryId = simpleMCPrimary.m_id;
                matchingDetails.m_nMatchedHits = simpleMatchedPfo.m_nMatchedHitsTotal;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::PrintMatchingOutput(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const
{
    std::cout << "---PROCESSED-MATCHING-OUTPUT--------------------------------------------------------------------" << std::endl;
    std::cout << "MinPrimaryGoodHits " << m_matchingMinPrimaryHits << ", MinHitsForGoodView " << m_matchingMinHitsForGoodView << ", MinPrimaryGoodViews " << m_matchingMinPrimaryGoodViews << std::endl;
    std::cout << "UseSmallPrimaries " << m_useSmallPrimaries << ", MinSharedHits " << m_matchingMinSharedHits << ", MinCompleteness " << m_matchingMinCompleteness << ", MinPurity " << m_matchingMinPurity << std::endl;

    bool isCorrect(true), isCalculable(false);

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
        const bool hasMatch(this->HasMatch(simpleMCPrimary, mapValue.second, matchingDetailsMap));
        const bool isTargetPrimary(this->IsGoodMCPrimary(simpleMCPrimary) && (NEUTRON != simpleMCPrimary.m_pdgCode));

        if (!hasMatch && !isTargetPrimary)
            continue;

        std::cout << std::endl << (!isTargetPrimary ? "(Non target) " : "")
                  << "Primary " << simpleMCPrimary.m_id << ", PDG " << simpleMCPrimary.m_pdgCode << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
                  << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << "),"
                  << " [nGood " << simpleMCPrimary.m_nGoodMCHitsTotal << " (" << simpleMCPrimary.m_nGoodMCHitsU << ", " << simpleMCPrimary.m_nGoodMCHitsV
                  << ", " << simpleMCPrimary.m_nGoodMCHitsW << ")]" << std::endl;

        if (NEUTRON != simpleMCPrimary.m_pdgCode)
            isCalculable = true;

        unsigned int nMatches(0);

        for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
        {
            if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
            {
                const bool isGoodMatch(this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo));

                if (isGoodMatch) ++nMatches;
                std::cout << "-" << (!isGoodMatch ? "(Below threshold) " : "") << "MatchedPfo " << simpleMatchedPfo.m_id;

                if (simpleMatchedPfo.m_parentId >= 0) std::cout << ", ParentPfo " << simpleMatchedPfo.m_parentId;

                std::cout << ", PDG " << simpleMatchedPfo.m_pdgCode << ", nMatchedHits " << simpleMatchedPfo.m_nMatchedHitsTotal
                          << " (" << simpleMatchedPfo.m_nMatchedHitsU << ", " << simpleMatchedPfo.m_nMatchedHitsV << ", " << simpleMatchedPfo.m_nMatchedHitsW << ")"
                          << ", nPfoHits " << simpleMatchedPfo.m_nPfoHitsTotal
                          << " (" << simpleMatchedPfo.m_nPfoHitsU << ", " << simpleMatchedPfo.m_nPfoHitsV << ", " << simpleMatchedPfo.m_nPfoHitsW << ")" << std::endl;
            }
        }

        if (isTargetPrimary && (1 != nMatches))
            isCorrect = false;
    }

    std::cout << std::endl << "Is correct? " << (isCorrect && isCalculable) << std::endl;
    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MONITORING
void EventValidationAlgorithm::VisualizeMatchingOutput(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector,
    const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const
{
    std::cout << "---VISUALIZE-MATCHING-OUTPUT--------------------------------------------------------------------" << std::endl;
    std::cout << "MinPrimaryGoodHits " << m_matchingMinPrimaryHits << ", MinHitsForGoodView " << m_matchingMinHitsForGoodView << ", MinPrimaryGoodViews " << m_matchingMinPrimaryGoodViews << std::endl;
    std::cout << "UseSmallPrimaries " << m_useSmallPrimaries << ", MinSharedHits " << m_matchingMinSharedHits << ", MinCompleteness " << m_matchingMinCompleteness << ", MinPurity " << m_matchingMinPurity << std::endl;

    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), m_visualizeGaps, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
    const HitTypeVector hitTypeVector{TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

    for (const HitType hitType : hitTypeVector)
    {
        const std::string hitTypeString((TPC_VIEW_U == hitType) ? "_U" : (TPC_VIEW_V == hitType) ? "_V" : "_W");

        if (m_visualizeVertices)
            this->VisualizeVertexMatches(mcNeutrinoVector, recoNeutrinoVector, hitType);

        for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
        {
            const SimpleMCPrimary &simpleMCPrimary(mapValue.first);
            bool hasMatch(this->HasMatch(simpleMCPrimary, mapValue.second, matchingDetailsMap));
            const bool isTargetPrimary(this->IsGoodMCPrimary(simpleMCPrimary) && (NEUTRON != simpleMCPrimary.m_pdgCode));

            if (!hasMatch && !isTargetPrimary)
                continue;

            PfoSet belowThresholdPfos;
            PfoVector primaryMatchedPfos;

            for (const SimpleMatchedPfo &simpleMatchedPfo : mapValue.second)
            {
                if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
                {
                    primaryMatchedPfos.push_back(simpleMatchedPfo.m_pPandoraAddress);

                    if (!this->IsGoodMatch(simpleMCPrimary, simpleMatchedPfo))
                        belowThresholdPfos.insert(simpleMatchedPfo.m_pPandoraAddress);
                }
            }

            std::string name("UNKNOWN_"); Color color(BLACK);
            this->GetPrimaryDetails(simpleMCPrimary, mcPrimaryMatchingMap, name, color);
            name.insert(0, !isTargetPrimary ? "NonTarget" : "");
            name += hitTypeString;

            if (isTargetPrimary && primaryMatchedPfos.empty())
            {
                const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), simpleMCPrimary.m_vertex, hitType));
                const CartesianVector endpointPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), simpleMCPrimary.m_endpoint, hitType));
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertexPosition2D, "MissingPrimaryVtx_" + name, RED, 1);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpointPosition2D, "MissingPrimaryEnd_" + name, RED, 1);
            }

            for (const ParticleFlowObject *const pPrimaryPfo : primaryMatchedPfos)
            {
                const bool isSplit((static_cast<int>(primaryMatchedPfos.size()) - static_cast<int>(belowThresholdPfos.size())) > 1);
                const bool isBestMatch(pPrimaryPfo == primaryMatchedPfos.front());
                const bool isBelowThreshold(belowThresholdPfos.count(pPrimaryPfo));
                const std::string prefix(isBelowThreshold ? "BelowThreshold" : !isSplit ? "Matched" : isBestMatch ? "BestFragment" : "Fragment");

                PfoList allPfos;
                LArPfoHelper::GetAllDownstreamPfos(pPrimaryPfo, allPfos);
                PfoVector allSortedPfos(allPfos.begin(), allPfos.end());
                std::sort(allSortedPfos.begin(), allSortedPfos.end(), LArPfoHelper::SortByNHits);

                for (const ParticleFlowObject *const pPfo : allSortedPfos)
                {
                    ClusterList clusterList;
                    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);
                    const std::string hierarchy((pPfo != pPrimaryPfo) ? "Daughter_" : "");
                    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &clusterList, prefix + "Clusters_" + hierarchy + name, color);

                    if (!isBelowThreshold && isSplit && !isBestMatch && !clusterList.empty())
                    {
                        CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
                        LArClusterHelper::GetExtremalCoordinates(*(clusterList.begin()), innerCoordinate, outerCoordinate);
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerCoordinate, "SplitMarker_" + hierarchy + name, MAGENTA, 1);
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerCoordinate, "SplitMarker_" + hierarchy + name, MAGENTA, 1);
                    }
                }
            }
        }

        if (m_visualizeRemnants)
            this->VisualizeRemnants(hitType);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::VisualizeVertexMatches(const MCParticleVector &mcNeutrinoVector, const PfoVector &recoNeutrinoVector, const HitType hitType) const
{
    const Vertex *pBestVertex(nullptr);
    float closestDistance(m_vertexVisualizationDeltaR);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinoVector)
        {
            const Vertex *const pVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
            const float distance((pMCNeutrino->GetEndpoint() - pVertex->GetPosition()).GetMagnitude());

            if (distance < closestDistance)
            {
                closestDistance = distance;
                pBestVertex = pVertex;
            }
        }
    }

    const std::string hitTypeString((TPC_VIEW_U == hitType) ? "_U" : (TPC_VIEW_V == hitType) ? "_V" : "_W");
    const bool isGood((1 == mcNeutrinoVector.size()) && (1 == recoNeutrinoVector.size()) && pBestVertex);
    const Color color(isGood ? GRAY : ORANGE);

    unsigned int displayIndex(0);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        const CartesianVector mcVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCNeutrino->GetEndpoint(), hitType));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertex2D, "MCNeutrinoVertex_" + TypeToString(displayIndex++) + hitTypeString, color, 1);
    }

    displayIndex = 0;

    for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinoVector)
    {
        const Vertex *const pThisVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
        const bool isThisGood(isGood || ((1 == mcNeutrinoVector.size()) && (pThisVertex == pBestVertex)));
        const std::string thisPrefix(isThisGood ? "Good" : "Displaced");
        const Color thisColor(isThisGood ? CYAN : DARKVIOLET);

        const CartesianVector recoVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pThisVertex->GetPosition(), hitType));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &recoVertex2D, thisPrefix + "RecoNeutrinoVertex_" + TypeToString(displayIndex++) + hitTypeString, thisColor, 1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetPrimaryDetails(const SimpleMCPrimary &thisSimpleMCPrimary, const MCPrimaryMatchingMap &mcPrimaryMatchingMap,
    std::string &name, Color &color) const
{
    // ATTN: Relies on fact that mcPrimaryMatchingMap is sorted by number of true hits
    unsigned int nLikeParticles(0);

    for (const MCPrimaryMatchingMap::value_type &mapValue : mcPrimaryMatchingMap)
    {
        const SimpleMCPrimary &mapSimpleMCPrimary(mapValue.first);

        if (thisSimpleMCPrimary.m_pdgCode == mapSimpleMCPrimary.m_pdgCode)
            ++nLikeParticles;

        if (thisSimpleMCPrimary.m_id == mapSimpleMCPrimary.m_id)
        {
            const std::string index(TypeToString(nLikeParticles));

            switch (mapSimpleMCPrimary.m_pdgCode)
            {
                case MU_MINUS: {name = "MU_MINUS" + index; color = RED;     return;}
                case MU_PLUS:  {name = "MU_PLUS"  + index; color = RED;     return;}
                case E_MINUS:  {name = "E_MINUS"  + index; color = ORANGE;  return;}
                case E_PLUS:   {name = "E_PLUS"   + index; color = ORANGE;  return;}
                case PROTON:   {name = "PROTON"   + index; color = (nLikeParticles % 2 == 0) ? VIOLET : BLUE;  return;}
                case PI_MINUS: {name = "PI_MINUS" + index; color = MAGENTA; return;}
                case PI_PLUS:  {name = "PI_PLUS"  + index; color = MAGENTA; return;}
                case PHOTON:   {name = "PHOTON"   + index; color = (nLikeParticles % 2 == 0) ? TEAL   : GREEN; return;}
                case NEUTRON:  {name = "NEUTRON"  + index; color = CYAN;    return;}
                default: return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::VisualizeRemnants(const HitType hitType) const
{
    const std::string hitTypeString((TPC_VIEW_U == hitType) ? "_U" : (TPC_VIEW_V == hitType) ? "_V" : "_W");

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    CaloHitList availableHits;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (PandoraContentApi::IsAvailable(*this, pCaloHit) && (hitType == pCaloHit->GetHitType()))
            availableHits.push_back(pCaloHit);
    }

    if (!availableHits.empty())
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &availableHits, "RemnantHits" + hitTypeString, LIGHTCYAN);

    CaloHitList isolatedHits;
    ClusterList availableClusters;

    for (const std::string &clusterListName : m_clusterListNames)
    {
        const ClusterList *pClusterList = nullptr;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList)
            continue;

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (PandoraContentApi::IsAvailable(*this, pCluster) && (hitType == LArClusterHelper::GetClusterHitType(pCluster)))
                availableClusters.push_back(pCluster);

            if (!PandoraContentApi::IsAvailable(*this, pCluster) && (hitType == LArClusterHelper::GetClusterHitType(pCluster)))
                isolatedHits.insert(isolatedHits.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        }
    }

    if (!availableClusters.empty())
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &availableClusters, "RemnantClusters" + hitTypeString, GRAY);

    if (!isolatedHits.empty())
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &isolatedHits, "IsolatedHitsInParticles" + hitTypeString, LIGHTYELLOW);
}
#endif
//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsGoodMCPrimary(const SimpleMCPrimary &simpleMCPrimary) const
{
    if (simpleMCPrimary.m_nGoodMCHitsTotal < m_matchingMinPrimaryHits)
        return false;

    int nGoodViews(0);
    if (simpleMCPrimary.m_nGoodMCHitsU >= m_matchingMinHitsForGoodView) ++nGoodViews;
    if (simpleMCPrimary.m_nGoodMCHitsV >= m_matchingMinHitsForGoodView) ++nGoodViews;
    if (simpleMCPrimary.m_nGoodMCHitsW >= m_matchingMinHitsForGoodView) ++nGoodViews;

    if (nGoodViews < m_matchingMinPrimaryGoodViews)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::HasMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfoList &simpleMatchedPfoList,
    const MatchingDetailsMap &matchingDetailsMap) const
{
    for (const SimpleMatchedPfo &simpleMatchedPfo : simpleMatchedPfoList)
    {
        if (matchingDetailsMap.count(simpleMatchedPfo.m_id) && (simpleMCPrimary.m_id == matchingDetailsMap.at(simpleMatchedPfo.m_id).m_matchedPrimaryId))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsGoodMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfo &simpleMatchedPfo) const
{
    const float purity((simpleMatchedPfo.m_nPfoHitsTotal > 0) ? static_cast<float>(simpleMatchedPfo.m_nMatchedHitsTotal) / static_cast<float>(simpleMatchedPfo.m_nPfoHitsTotal) : 0.f);
    const float completeness((simpleMCPrimary.m_nMCHitsTotal > 0) ? static_cast<float>(simpleMatchedPfo.m_nMatchedHitsTotal) / static_cast<float>(simpleMCPrimary.m_nMCHitsTotal) : 0.f);

    return ((simpleMatchedPfo.m_nMatchedHitsTotal >= m_matchingMinSharedHits) && (purity >= m_matchingMinPurity) && (completeness >= m_matchingMinCompleteness));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetPfoIdMap(const PfoList &pfoList, PfoIdMap &pfoIdMap) const
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
    if (lhs.m_nGoodMCHitsTotal != rhs.m_nGoodMCHitsTotal)
        return (lhs.m_nGoodMCHitsTotal > rhs.m_nGoodMCHitsTotal);

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

bool EventValidationAlgorithm::SortRecoNeutrinos(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    if (!LArPfoHelper::IsNeutrino(pLhs) || !LArPfoHelper::IsNeutrino(pRhs))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PfoList downstreamPfosLhs, downstreamPfosRhs;
    LArPfoHelper::GetAllDownstreamPfos(pLhs, downstreamPfosLhs);
    LArPfoHelper::GetAllDownstreamPfos(pRhs, downstreamPfosRhs);

    // If just left with the neutrino pfos themselves
    if ((1 == downstreamPfosLhs.size()) && (1 == downstreamPfosRhs.size()))
    {
        // ATTN Not a good pair of tie-breakers, but this should be rare (ideally shouldn't have any neutrinos without daughter pfos)
        if (!pLhs->GetVertexList().empty() && !pRhs->GetVertexList().empty())
            return ((*(pLhs->GetVertexList().begin()))->GetPosition().GetZ() < (*(pRhs->GetVertexList().begin()))->GetPosition().GetZ());

        return (pLhs->GetParticleId() < pRhs->GetParticleId());
    }

    PfoVector pfoVectorLhs(downstreamPfosLhs.begin(), downstreamPfosLhs.end());
    PfoVector pfoVectorRhs(downstreamPfosRhs.begin(), downstreamPfosRhs.end());
    std::sort(pfoVectorLhs.begin(), pfoVectorLhs.end(), LArPfoHelper::SortByNHits);
    std::sort(pfoVectorRhs.begin(), pfoVectorRhs.end(), LArPfoHelper::SortByNHits);

    if (pfoVectorLhs.empty() || pfoVectorRhs.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return (LArPfoHelper::SortByNHits(pfoVectorLhs.front(), pfoVectorRhs.front()));
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
    m_nGoodMCHitsTotal(0),
    m_nGoodMCHitsU(0),
    m_nGoodMCHitsV(0),
    m_nGoodMCHitsW(0),
    m_energy(0.f),
    m_momentum(0.f, 0.f, 0.f),
    m_vertex(-1.f, -1.f, -1.f),
    m_endpoint(-1.f, -1.f, -1.f),
    m_nMatchedPfos(0),
    m_pPandoraAddress(nullptr)
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
    m_pPandoraAddress(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::MatchingDetails::MatchingDetails() :
    m_matchedPrimaryId(-1),
    m_nMatchedHits(0)
{
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
        "IntegrateOverSlices", m_integrateOverSlices));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoInducedOnly", m_neutrinoInducedOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrimaryPfosOnly", m_primaryPfosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CollapseToPrimaryPfos", m_collapseToPrimaryPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeMatching", m_visualizeMatching));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeVertices", m_visualizeVertices));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
         "VisualizeRemnants", m_visualizeRemnants));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeGaps", m_visualizeGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPrimaryHits", m_matchingMinPrimaryHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinHitsForGoodView", m_matchingMinHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPrimaryGoodViews", m_matchingMinPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseSmallPrimaries", m_useSmallPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinSharedHits", m_matchingMinSharedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinCompleteness", m_matchingMinCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPurity", m_matchingMinPurity));

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
