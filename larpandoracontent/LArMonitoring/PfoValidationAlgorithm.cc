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

PfoValidationAlgorithm::PfoValidationAlgorithm()
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

    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamParticle, beamMCParticlesToGoodHitsMap);
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsCosmicRay, crMCParticlesToGoodHitsMap);

    LArMCParticleHelper::MCContributionMapVector mcParticlesToGoodHitsMaps;
    mcParticlesToGoodHitsMaps.push_back(nuMCParticlesToGoodHitsMap);
    mcParticlesToGoodHitsMaps.push_back(beamMCParticlesToGoodHitsMap);
    mcParticlesToGoodHitsMaps.push_back(crMCParticlesToGoodHitsMap);

    // Get the mappings detailing the hits shared between Pfos and reconstructable MCParticles
    pandora::PfoList finalStatePfos;

    // TODO use helper function in LArMonitoringHelper, not sure if it currently give the right output
    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }

    LArMCParticleHelper::PfoContributionMap pfoToReconstructable2DHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, mcParticlesToGoodHitsMaps, pfoToReconstructable2DHitsMap);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToReconstructable2DHitsMap, mcParticlesToGoodHitsMaps, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

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

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodHits", m_parameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForGoodView", m_parameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodViews", m_parameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_parameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_parameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_parameters.m_minHitSharingFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
