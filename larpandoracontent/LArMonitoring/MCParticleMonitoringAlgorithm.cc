/**
 *  @file   larpandoracontent/LArMonitoring/MCParticleMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the mc particle monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/MCParticleMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{

MCParticleMonitoringAlgorithm::MCParticleMonitoringAlgorithm() :
    m_neutrinoInducedOnly(false),
    m_minHitsForDisplay(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCParticleMonitoringAlgorithm::Run()
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // MC primary records, with details of all hits associated with all particles in decay chain
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);

    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    LArMonitoringHelper::CaloHitToMCMap hitToPrimaryMCMap;
    LArMonitoringHelper::MCContributionMap mcPrimaryToTrueHitListMap;
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcPrimaryToTrueHitListMap);

    SimpleMCParticleList simpleMCPrimaryList;
    this->GetSimpleMCParticleList(mcPrimaryVector, mcPrimaryToTrueHitListMap, simpleMCPrimaryList);

    // MC particle records, with details of hits associated to individual mc particles
    const MCParticleVector mcParticleVector(pMCParticleList->begin(), pMCParticleList->end());

    LArMCParticleHelper::MCRelationMap mcParticleToSelfMap;
    for (const MCParticle *const pMCParticle : mcParticleVector) mcParticleToSelfMap[pMCParticle] = pMCParticle;

    LArMonitoringHelper::CaloHitToMCMap hitToMCMap;
    LArMonitoringHelper::MCContributionMap mcToTrueHitListMap;
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcParticleToSelfMap, hitToMCMap, mcToTrueHitListMap);    

    SimpleMCParticleList simpleMCParticleList;
    this->GetSimpleMCParticleList(mcParticleVector, mcToTrueHitListMap, simpleMCParticleList);

    SimpleMCParticleMap simpleMCParticleMap;
    for (const SimpleMCParticle &simpleMCParticle : simpleMCParticleList) simpleMCParticleMap[simpleMCParticle.m_pPandoraAddress] = simpleMCParticle;

    // Display full hierarchy for each mc neutrino
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoVector);
    this->PrintAllOutput(mcNeutrinoVector, simpleMCPrimaryList, simpleMCParticleMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCParticleMonitoringAlgorithm::GetSimpleMCParticleList(const MCParticleVector &mcPrimaryVector, const LArMonitoringHelper::MCContributionMap &mcToTrueHitListMap,
    SimpleMCParticleList &simpleMCParticleList) const
{
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (m_neutrinoInducedOnly && !LArMCParticleHelper::IsNeutrinoInduced(pMCPrimary))
            continue;

        SimpleMCParticle simpleMCParticle;
        // ATTN simpleMCParticle.m_id assigned later, after sorting
        simpleMCParticle.m_pPandoraAddress = pMCPrimary;
        simpleMCParticle.m_pdgCode = pMCPrimary->GetParticleId();
        simpleMCParticle.m_energy = pMCPrimary->GetEnergy();
        simpleMCParticle.m_momentum = pMCPrimary->GetMomentum();
        simpleMCParticle.m_vertex = pMCPrimary->GetVertex();
        simpleMCParticle.m_endpoint = pMCPrimary->GetEndpoint();

        LArMonitoringHelper::MCContributionMap::const_iterator trueHitsIter = mcToTrueHitListMap.find(pMCPrimary);

        if (mcToTrueHitListMap.end() != trueHitsIter)
        {
            const CaloHitList &caloHitList(trueHitsIter->second);
            simpleMCParticle.m_nMCHitsTotal = caloHitList.size();
            simpleMCParticle.m_nMCHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList);
            simpleMCParticle.m_nMCHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList);
            simpleMCParticle.m_nMCHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList);
        }

        simpleMCParticleList.push_back(simpleMCParticle);
    }

    std::sort(simpleMCParticleList.begin(), simpleMCParticleList.end(), MCParticleMonitoringAlgorithm::SortSimpleMCParticles);

    int mcPrimaryId(0);
    for (SimpleMCParticle &simpleMCParticle : simpleMCParticleList)
        simpleMCParticle.m_id = mcPrimaryId++;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCParticleMonitoringAlgorithm::PrintAllOutput(const MCParticleVector &mcNeutrinoVector, const SimpleMCParticleList &simpleMCPrimaryList,
    const SimpleMCParticleMap &simpleMCParticleMap) const
{
    std::cout << "---MC-PARTICLE-MONITORING-----------------------------------------------------------------------" << std::endl;

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);
        std::cout << "---MCNeutrino, PDG " << pMCNeutrino->GetParticleId() << ", Energy " << pMCNeutrino->GetEnergy() << ", Nuance " << (pLArMCNeutrino ? pLArMCNeutrino->GetNuanceCode() : -1) << std::endl;

        for (const SimpleMCParticle &simpleMCPrimary : simpleMCPrimaryList)
        {
            const MCParticle *const pMCPrimary(simpleMCPrimary.m_pPandoraAddress);

            if (pMCNeutrino != LArMCParticleHelper::GetParentNeutrino(pMCPrimary))
                continue;

            std::cout << std::endl << "--Primary " << simpleMCPrimary.m_id << ", MCPDG " << simpleMCPrimary.m_pdgCode << ", Energy " << simpleMCPrimary.m_energy
                      << ", nMCHits " << simpleMCPrimary.m_nMCHitsTotal
                      << " (" << simpleMCPrimary.m_nMCHitsU << ", " << simpleMCPrimary.m_nMCHitsV << ", " << simpleMCPrimary.m_nMCHitsW << ")" << std::endl;

            this->PrintMCParticle(pMCPrimary, simpleMCParticleMap, 1);
        }

        std::cout << std::endl;
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCParticleMonitoringAlgorithm::PrintMCParticle(const MCParticle *const pMCParticle, const SimpleMCParticleMap &simpleMCParticleMap,
    const int depth) const
{
    const SimpleMCParticle &simpleMCParticle(simpleMCParticleMap.at(pMCParticle));

    if (simpleMCParticle.m_nMCHitsTotal >= m_minHitsForDisplay)
    {
        if (depth > 1)
        {
            for (int iDepth = 1; iDepth < depth - 1; ++iDepth) std::cout << "   ";
            std::cout << "\\_ ";
        }

        std::cout << "MCPDG " << simpleMCParticle.m_pdgCode << ", Energy " << simpleMCParticle.m_energy << ", nMCHits " << simpleMCParticle.m_nMCHitsTotal
                  << " (" << simpleMCParticle.m_nMCHitsU << ", " << simpleMCParticle.m_nMCHitsV << ", " << simpleMCParticle.m_nMCHitsW << ")" << std::endl;
    }

    for (const MCParticle *const pDaughterParticle : pMCParticle->GetDaughterList())
        this->PrintMCParticle(pDaughterParticle, simpleMCParticleMap, depth + 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MCParticleMonitoringAlgorithm::SortSimpleMCParticles(const SimpleMCParticle &lhs, const SimpleMCParticle &rhs)
{
    if (lhs.m_nMCHitsTotal != rhs.m_nMCHitsTotal)
        return (lhs.m_nMCHitsTotal > rhs.m_nMCHitsTotal);

    return (lhs.m_energy > rhs.m_energy);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleMonitoringAlgorithm::SimpleMCParticle::SimpleMCParticle() :
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
    m_pPandoraAddress(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MCParticleMonitoringAlgorithm::SimpleMCParticle::operator<(const SimpleMCParticle &rhs) const
{
    if (this == &rhs)
        return false;

    return (m_id < rhs.m_id);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoInducedOnly", m_neutrinoInducedOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForDisplay", m_minHitsForDisplay));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
