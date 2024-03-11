/**
 *  @file   larpandoracontent/LArMonitoring/PfoValidationAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/PfoValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoValidationAlgorithm::PfoValidationAlgorithm() :
    m_nMatchesToShow(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoValidationAlgorithm::Run()
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    // Identify reconstructable MCParticles, and get mappings to their good hits
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::MCContributionMap beamMCParticlesToGoodHitsMap;
    LArMCParticleHelper::MCContributionMap crMCParticlesToGoodHitsMap;

    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamParticle, beamMCParticlesToGoodHitsMap);
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsCosmicRay, crMCParticlesToGoodHitsMap);

    const LArMCParticleHelper::MCContributionMapVector mcParticlesToGoodHitsMaps{
        nuMCParticlesToGoodHitsMap, beamMCParticlesToGoodHitsMap, crMCParticlesToGoodHitsMap};

    // Get the mappings detailing the hits shared between Pfos and reconstructable MCParticles
    PfoList finalStatePfos;

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }

    LArMCParticleHelper::PfoContributionMap pfoToReconstructable2DHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(
        finalStatePfos, mcParticlesToGoodHitsMaps, pfoToReconstructable2DHitsMap, m_parameters.m_foldBackHierarchy);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(
        pfoToReconstructable2DHitsMap, mcParticlesToGoodHitsMaps, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Print the monte-carlo information for this event
    MCParticleVector orderedMCParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector(mcParticlesToGoodHitsMaps, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("MC : Reconstructable neutrino final state particles");
    LArMonitoringHelper::PrintMCParticleTable(nuMCParticlesToGoodHitsMap, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("MC : Reconstructable primary beam particles");
    LArMonitoringHelper::PrintMCParticleTable(beamMCParticlesToGoodHitsMap, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("MC : Reconstructable primary cosmic-rays");
    LArMonitoringHelper::PrintMCParticleTable(crMCParticlesToGoodHitsMap, orderedMCParticleVector);

    // Print the pfo information for this event
    PfoVector orderedPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToReconstructable2DHitsMap, orderedPfoVector);

    LArFormattingHelper::PrintHeader("Reco : Primary Pfos");
    LArMonitoringHelper::PrintPfoTable(pfoToReconstructable2DHitsMap, orderedPfoVector);

    // Print the raw matching between Pfos and MCParticles
    LArFormattingHelper::PrintHeader("Raw Reco vs. MC matching");
    LArMonitoringHelper::PrintMatchingTable(orderedPfoVector, orderedMCParticleVector, mcParticleToPfoHitSharingMap, m_nMatchesToShow);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodHits", m_parameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitsForGoodView", m_parameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodViews", m_parameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectInputHits", m_parameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxPhotonPropagation", m_parameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitSharingFraction", m_parameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NumMatchesToShow", m_nMatchesToShow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
