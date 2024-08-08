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
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/MCParticleMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{


MCParticleMonitoringAlgorithm::MCParticleMonitoringAlgorithm() : 
    m_useTrueNeutrinosOnly(false), 
    m_minHitsForDisplay(1), 
    m_minPrimaryGoodHits(3), 
    m_minHitsForGoodView(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCParticleMonitoringAlgorithm::Run()
{
    std::cout << "---MC-PARTICLE-MONITORING-----------------------------------------------------------------------" << std::endl;
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = m_minPrimaryGoodHits;
    parameters.m_minHitsForGoodView = m_minHitsForGoodView;
    parameters.m_minHitSharingFraction = 0.f;

    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::MCContributionMap beamMCParticlesToGoodHitsMap;
    LArMCParticleHelper::MCContributionMap crMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    if (!m_useTrueNeutrinosOnly)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, beamMCParticlesToGoodHitsMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, crMCParticlesToGoodHitsMap);
    }

    if (!nuMCParticlesToGoodHitsMap.empty())
    {
        std::cout << std::endl << "BeamNeutrinos: " << std::endl;
        this->PrintPrimaryMCParticles(nuMCParticlesToGoodHitsMap);
    }

    if (!beamMCParticlesToGoodHitsMap.empty())
    {
        std::cout << std::endl << "BeamParticles: " << std::endl;
        this->PrintPrimaryMCParticles(beamMCParticlesToGoodHitsMap);
    }

    if (!crMCParticlesToGoodHitsMap.empty())
    {
        std::cout << std::endl << "CosmicRays: " << std::endl;
        this->PrintPrimaryMCParticles(crMCParticlesToGoodHitsMap);
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCParticleMonitoringAlgorithm::PrintPrimaryMCParticles(const LArMCParticleHelper::MCContributionMap &mcContributionMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcContributionMap}, mcPrimaryVector);

    unsigned int index(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const CaloHitList &caloHitList(mcContributionMap.at(pMCPrimary));
        const LArMCParticle *pLArMC{dynamic_cast<const LArMCParticle *>(pMCPrimary)};
        if (caloHitList.size() >= m_minHitsForDisplay)
        {
            std::cout << std::endl
                      << "--Primary " << index << ", MCPDG " << pMCPrimary->GetParticleId() << ", Energy " << pMCPrimary->GetEnergy()
                      << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude() << ", nMCHits "
                      << caloHitList.size() << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList) << ", "
                      << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList) << ", "
                      << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList) << " , Process: " << pLArMC->GetProcess() << ")" << std::endl;

	    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
            LArMCParticleHelper::CaloHitToMCMap caloHitToPrimaryMCMap;
            LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
            LArMCParticleHelper::GetMCParticleToCaloHitMatches(&caloHitList, mcToPrimaryMCMap, caloHitToPrimaryMCMap, mcToTrueHitListMap);
            this->PrintMCParticle(pMCPrimary, mcToTrueHitListMap, 1);
        }

        ++index;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCParticleMonitoringAlgorithm::PrintMCParticle(
    const MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap, const int depth) const
{
    const CaloHitList &caloHitList(mcToTrueHitListMap.count(pMCParticle) ? mcToTrueHitListMap.at(pMCParticle) : CaloHitList());

    if (caloHitList.size() >= m_minHitsForDisplay)
    {
        if (depth > 1)
        {
            for (int iDepth = 1; iDepth < depth - 1; ++iDepth)
                std::cout << "   ";
            std::cout << "\\_ ";
        }

        std::cout << "MCPDG " << pMCParticle->GetParticleId() << ", Energy " << pMCParticle->GetEnergy() << ", Dist. "
                  << (pMCParticle->GetEndpoint() - pMCParticle->GetVertex()).GetMagnitude() << ", nMCHits " << caloHitList.size() << " ("
                  << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList) << ", "
                  << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList) << ", "
                  << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList) << ")" << std::endl;
    }

    for (const MCParticle *const pDaughterParticle : pMCParticle->GetDaughterList())
        this->PrintMCParticle(pDaughterParticle, mcToTrueHitListMap, depth + 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitsForDisplay", m_minHitsForDisplay));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodHits", m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitsForGoodView", m_minHitsForGoodView));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
