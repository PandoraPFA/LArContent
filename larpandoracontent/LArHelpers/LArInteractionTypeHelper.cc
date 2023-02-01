/**
 *  @file   larpandoracontent/LArHelpers/LArInteractionTypeHelper.cc
 *
 *  @brief  Implementation of the interaction type helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

using namespace pandora;

namespace lar_content
{

LArInteractionTypeHelper::InteractionType LArInteractionTypeHelper::GetInteractionType(const MCParticleList &mcPrimaryList)
{
    if (mcPrimaryList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    InteractionParameters parameters;
    LArInteractionTypeHelper::SetInteractionParameters(mcPrimaryList, parameters);

    const InteractionType cosmicRayHypothesis(LArInteractionTypeHelper::CosmicRayHypothesis(mcPrimaryList, parameters));

    if (OTHER_INTERACTION != cosmicRayHypothesis)
        return cosmicRayHypothesis;

    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsBeamParticle(mcPrimaryList.front()))
    {
        if (1 == parameters.m_nMuons)
            return BEAM_PARTICLE_MU;
        if (1 == parameters.m_nProtons)
            return BEAM_PARTICLE_P;
        if (1 == parameters.m_nElectrons)
            return BEAM_PARTICLE_E;
        if (1 == parameters.m_nPhotons)
            return BEAM_PARTICLE_PHOTON;
        if (1 == parameters.m_nPiPlus)
            return BEAM_PARTICLE_PI_PLUS;
        if (1 == parameters.m_nPiMinus)
            return BEAM_PARTICLE_PI_MINUS;
        if (1 == parameters.m_nKaonPlus)
            return BEAM_PARTICLE_KAON_PLUS;
        if (1 == parameters.m_nKaonMinus)
            return BEAM_PARTICLE_KAON_MINUS;
        else
            return BEAM_PARTICLE_OTHER;
    }

    const MCParticle *pMCNeutrino(nullptr);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (!LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) ||
            (pMCNeutrino && (pMCNeutrino != LArMCParticleHelper::GetParentMCParticle(pMCPrimary))))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pMCNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCPrimary);
    }

    if (!pMCNeutrino)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const int nuNuanceCode(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));

    if (1001 == nuNuanceCode)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons))
            return CCQEL_MU;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons))
            return CCQEL_MU_P;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons))
            return CCQEL_MU_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons))
            return CCQEL_MU_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons))
            return CCQEL_MU_P_P_P_P;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons))
            return CCQEL_MU_P_P_P_P_P;

        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (0 == parameters.m_nProtons))
            return CCQEL_E;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (1 == parameters.m_nProtons))
            return CCQEL_E_P;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (2 == parameters.m_nProtons))
            return CCQEL_E_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (3 == parameters.m_nProtons))
            return CCQEL_E_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (4 == parameters.m_nProtons))
            return CCQEL_E_P_P_P_P;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (5 == parameters.m_nProtons))
            return CCQEL_E_P_P_P_P_P;
    }

    if (1002 == nuNuanceCode)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons))
            return NCQEL_P;
        if ((2 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons))
            return NCQEL_P_P;
        if ((3 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons))
            return NCQEL_P_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons))
            return NCQEL_P_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons))
            return NCQEL_P_P_P_P_P;
    }

    if ((nuNuanceCode >= 1003) && (nuNuanceCode <= 1005))
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons))
            return CCRES_MU;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons))
            return CCRES_MU_P;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons))
            return CCRES_MU_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons))
            return CCRES_MU_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons))
            return CCRES_MU_P_P_P_P;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons))
            return CCRES_MU_P_P_P_P_P;

        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_MU_PIPLUS;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_MU_P_PIPLUS;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_MU_P_P_PIPLUS;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_MU_P_P_P_PIPLUS;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_MU_P_P_P_P_PIPLUS;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_MU_P_P_P_P_P_PIPLUS;

        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_MU_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_MU_P_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_MU_P_P_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_MU_P_P_P_PHOTON;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_MU_P_P_P_P_PHOTON;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_MU_P_P_P_P_P_PHOTON;

        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_MU_PIZERO;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_MU_P_PIZERO;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_MU_P_P_PIZERO;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_MU_P_P_P_PIZERO;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_MU_P_P_P_P_PIZERO;
        if ((8 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_MU_P_P_P_P_P_PIZERO;

        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (0 == parameters.m_nProtons))
            return CCRES_E;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (1 == parameters.m_nProtons))
            return CCRES_E_P;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (2 == parameters.m_nProtons))
            return CCRES_E_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (3 == parameters.m_nProtons))
            return CCRES_E_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (4 == parameters.m_nProtons))
            return CCRES_E_P_P_P_P;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (5 == parameters.m_nProtons))
            return CCRES_E_P_P_P_P_P;

        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_E_PIPLUS;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_E_P_PIPLUS;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_E_P_P_PIPLUS;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_E_P_P_P_PIPLUS;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_E_P_P_P_P_PIPLUS;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCRES_E_P_P_P_P_P_PIPLUS;

        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_E_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_E_P_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_E_P_P_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_E_P_P_P_PHOTON;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_E_P_P_P_P_PHOTON;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCRES_E_P_P_P_P_P_PHOTON;

        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (0 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_E_PIZERO;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (1 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_E_P_PIZERO;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (2 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_E_P_P_PIZERO;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (3 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_E_P_P_P_PIZERO;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (4 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_E_P_P_P_P_PIZERO;
        if ((8 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons) && (5 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCRES_E_P_P_P_P_P_PIZERO;
    }

    if ((nuNuanceCode >= 1006) && (nuNuanceCode <= 1009))
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons))
            return NCRES_P;
        if ((2 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons))
            return NCRES_P_P;
        if ((3 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons))
            return NCRES_P_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons))
            return NCRES_P_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons))
            return NCRES_P_P_P_P_P;

        if ((1 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCRES_PIPLUS;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCRES_P_PIPLUS;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCRES_P_P_PIPLUS;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCRES_P_P_P_PIPLUS;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCRES_P_P_P_P_PIPLUS;
        if ((6 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCRES_P_P_P_P_P_PIPLUS;

        if ((1 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCRES_PIMINUS;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCRES_P_PIMINUS;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCRES_P_P_PIMINUS;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCRES_P_P_P_PIMINUS;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCRES_P_P_P_P_PIMINUS;
        if ((6 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCRES_P_P_P_P_P_PIMINUS;

        if ((1 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCRES_PHOTON;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCRES_P_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCRES_P_P_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCRES_P_P_P_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCRES_P_P_P_P_PHOTON;
        if ((6 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCRES_P_P_P_P_P_PHOTON;

        if ((2 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCRES_PIZERO;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCRES_P_PIZERO;
        if ((4 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCRES_P_P_PIZERO;
        if ((5 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCRES_P_P_P_PIZERO;
        if ((6 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCRES_P_P_P_P_PIZERO;
        if ((7 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCRES_P_P_P_P_P_PIZERO;
    }

    if (1091 == nuNuanceCode)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons))
            return CCDIS_MU;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons))
            return CCDIS_MU_P;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons))
            return CCDIS_MU_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons))
            return CCDIS_MU_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons))
            return CCDIS_MU_P_P_P_P;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons))
            return CCDIS_MU_P_P_P_P_P;

        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCDIS_MU_PIPLUS;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCDIS_MU_P_PIPLUS;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCDIS_MU_P_P_PIPLUS;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCDIS_MU_P_P_P_PIPLUS;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCDIS_MU_P_P_P_P_PIPLUS;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return CCDIS_MU_P_P_P_P_P_PIPLUS;

        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCDIS_MU_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCDIS_MU_P_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_P_PHOTON;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_P_P_PHOTON;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_P_P_P_PHOTON;

        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (0 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCDIS_MU_PIZERO;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (1 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCDIS_MU_P_PIZERO;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (2 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_PIZERO;
        if ((6 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (3 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_P_PIZERO;
        if ((7 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (4 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_P_P_PIZERO;
        if ((8 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons) && (5 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return CCDIS_MU_P_P_P_P_P_PIZERO;
    }

    if (1092 == nuNuanceCode)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons))
            return NCDIS_P;
        if ((2 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons))
            return NCDIS_P_P;
        if ((3 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons))
            return NCDIS_P_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons))
            return NCDIS_P_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons))
            return NCDIS_P_P_P_P_P;

        if ((1 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCDIS_PIPLUS;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCDIS_P_PIPLUS;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCDIS_P_P_PIPLUS;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCDIS_P_P_P_PIPLUS;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCDIS_P_P_P_P_PIPLUS;
        if ((6 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiPlus))
            return NCDIS_P_P_P_P_P_PIPLUS;

        if ((1 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCDIS_PIMINUS;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCDIS_P_PIMINUS;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCDIS_P_P_PIMINUS;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCDIS_P_P_P_PIMINUS;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCDIS_P_P_P_P_PIMINUS;
        if ((6 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPiMinus))
            return NCDIS_P_P_P_P_P_PIMINUS;

        if ((1 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCDIS_PHOTON;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCDIS_P_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCDIS_P_P_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCDIS_P_P_P_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCDIS_P_P_P_P_PHOTON;
        if ((6 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return NCDIS_P_P_P_P_P_PHOTON;

        if ((2 == parameters.m_nNonNeutrons) && (0 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCDIS_PIZERO;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCDIS_P_PIZERO;
        if ((4 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCDIS_P_P_PIZERO;
        if ((5 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCDIS_P_P_P_PIZERO;
        if ((6 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCDIS_P_P_P_P_PIZERO;
        if ((7 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return NCDIS_P_P_P_P_P_PIZERO;
    }

    if (1096 == nuNuanceCode)
        return NCCOH;
    if (1097 == nuNuanceCode)
        return CCCOH;

    return OTHER_INTERACTION;
}

//------------------------------------------------------------------------------------------------------------------------------------------

InteractionDescriptor LArInteractionTypeHelper::GetInteractionDescriptor(const MCParticleList &mcPrimaryList)
{
    if (mcPrimaryList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    InteractionParameters parameters;
    LArInteractionTypeHelper::SetInteractionParameters(mcPrimaryList, parameters);

    const MCParticle *pMCNeutrino(nullptr);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (!LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) ||
            (pMCNeutrino && (pMCNeutrino != LArMCParticleHelper::GetParentMCParticle(pMCPrimary))))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pMCNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCPrimary);
    }

    if (!pMCNeutrino)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const int nuNuanceCode(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));

    if (1001 == nuNuanceCode)
    {
        return InteractionDescriptor(true, true, false, false, false, parameters.m_nMuons == 1, parameters.m_nElectrons == 1, 0, 0, 0,
            parameters.m_nProtons);
    }

    if (1002 == nuNuanceCode)
    {
        return InteractionDescriptor(false, true, false, false, false, false, false, 0, 0, 0, parameters.m_nProtons);
    }

    if ((nuNuanceCode >= 1003) && (nuNuanceCode <= 1005))
    {
        return InteractionDescriptor(true, false, true, false, false, parameters.m_nMuons == 1, parameters.m_nElectrons == 1,
            parameters.m_nPiPlus, parameters.m_nPiMinus, parameters.m_nPhotons, parameters.m_nProtons);
    }

    if ((nuNuanceCode >= 1006) && (nuNuanceCode <= 1009))
    {
        return InteractionDescriptor(false, false, true, false, false, parameters.m_nMuons == 1, parameters.m_nElectrons == 1,
            parameters.m_nPiPlus, parameters.m_nPiMinus, parameters.m_nPhotons, parameters.m_nProtons);
    }

    if (1091 == nuNuanceCode)
    {
        return InteractionDescriptor(true, false, false, true, false, parameters.m_nMuons == 1, parameters.m_nElectrons == 1,
            parameters.m_nPiPlus, parameters.m_nPiMinus, parameters.m_nPhotons, parameters.m_nProtons);
    }

    if (1092 == nuNuanceCode)
    {
        return InteractionDescriptor(false, false, false, true, false, parameters.m_nMuons == 1, parameters.m_nElectrons == 1,
            parameters.m_nPiPlus, parameters.m_nPiMinus, parameters.m_nPhotons, parameters.m_nProtons);
    }

    if (1096 == nuNuanceCode)
    {
        return InteractionDescriptor(false, false, false, false, true, parameters.m_nMuons == 1, parameters.m_nElectrons == 1, 0, 0, 0, 0);
    }
    if (1097 == nuNuanceCode)
    {
        return InteractionDescriptor(true, false, false, false, true, parameters.m_nMuons == 1, parameters.m_nElectrons == 1, 0, 0, 0, 0);
    }

    return InteractionDescriptor(false, false, false, false, false, false, false, 0, 0, 0, 0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType LArInteractionTypeHelper::GetTestBeamHierarchyInteractionType(const MCParticleList &mcPrimaryList)
{
    if (mcPrimaryList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    InteractionParameters parameters;
    LArInteractionTypeHelper::SetInteractionParameters(mcPrimaryList, parameters);

    const InteractionType cosmicRayHypothesis(LArInteractionTypeHelper::CosmicRayHypothesis(mcPrimaryList, parameters));

    if (OTHER_INTERACTION != cosmicRayHypothesis)
        return cosmicRayHypothesis;

    int targetParentId(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (LArMCParticleHelper::GetParentMCParticle(pMCPrimary) != pMCPrimary)
            continue;

        if (13 == std::fabs(pMCPrimary->GetParticleId()))
            --parameters.m_nMuons;
        else if (11 == std::fabs(pMCPrimary->GetParticleId()))
            --parameters.m_nElectrons;
        else if (2212 == std::fabs(pMCPrimary->GetParticleId()))
            --parameters.m_nProtons;
        else if (22 == pMCPrimary->GetParticleId())
            --parameters.m_nPhotons;
        else if (211 == pMCPrimary->GetParticleId())
            --parameters.m_nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId())
            --parameters.m_nPiMinus;
        else if (321 == pMCPrimary->GetParticleId())
            --parameters.m_nKaonPlus;
        else if (-321 == pMCPrimary->GetParticleId())
            --parameters.m_nKaonMinus;

        --parameters.m_nNonNeutrons;
        targetParentId = pMCPrimary->GetParticleId();
    }

    MCParticleSet piZero, kaon0L;
    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (22 == pMCPrimary->GetParticleId() && pMCPrimary->GetParentList().size() == 1)
        {
            const MCParticle *const pParentMCParticle(pMCPrimary->GetParentList().front());
            if (pParentMCParticle->GetParticleId() == 111 && !piZero.count(pParentMCParticle))
                piZero.insert(pParentMCParticle);
        }
        else if ((111 == pMCPrimary->GetParticleId() || 211 == std::fabs(pMCPrimary->GetParticleId())) && pMCPrimary->GetParentList().size() == 1)
        {
            const MCParticle *const pParentMCParticle(pMCPrimary->GetParentList().front());
            if (pParentMCParticle->GetParticleId() == 130 && !kaon0L.count(pParentMCParticle))
                kaon0L.insert(pParentMCParticle);
        }
    }

    parameters.m_nPiZero = piZero.size();
    parameters.m_nKaon0L = kaon0L.size();

    if (211 == targetParentId)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nPiPlus))
            return BEAM_PARTICLE_PI_PLUS_PI_PLUS;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nPiPlus) && (1 == parameters.m_nPhotons))
            return BEAM_PARTICLE_PI_PLUS_PI_PLUS_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nPiPlus) && (1 == parameters.m_nPiZero))
            return BEAM_PARTICLE_PI_PLUS_PI_PLUS_PIZERO;
        return BEAM_PARTICLE_PI_PLUS_COMPLEX;
    }
    else if (-211 == targetParentId)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nPiMinus))
            return BEAM_PARTICLE_PI_MINUS_PI_MINUS;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nPiMinus) && (1 == parameters.m_nPhotons))
            return BEAM_PARTICLE_PI_MINUS_PI_MINUS_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nPiMinus) && (1 == parameters.m_nPiZero))
            return BEAM_PARTICLE_PI_MINUS_PI_MINUS_PIZERO;
        return BEAM_PARTICLE_PI_MINUS_COMPLEX;
    }
    else if (2212 == targetParentId)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons))
            return BEAM_PARTICLE_P_P;
        if ((2 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_PHOTON;
        if ((3 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_PHOTON_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (3 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (1 == parameters.m_nProtons) && (4 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON_PHOTON;

        if ((2 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons))
            return BEAM_PARTICLE_P_P_P;
        if ((3 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_P_PHOTON;
        if ((4 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_P_PHOTON_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (2 == parameters.m_nProtons) && (3 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_P_PHOTON_PHOTON_PHOTON;

        if ((3 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons))
            return BEAM_PARTICLE_P_P_P_P;
        if ((4 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_P_P_PHOTON;
        if ((5 == parameters.m_nNonNeutrons) && (3 == parameters.m_nProtons) && (2 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_P_P_PHOTON_PHOTON;

        if ((4 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons))
            return BEAM_PARTICLE_P_P_P_P_P;
        if ((5 == parameters.m_nNonNeutrons) && (4 == parameters.m_nProtons) && (1 == parameters.m_nPhotons))
            return BEAM_PARTICLE_P_P_P_P_P_PHOTON;

        if ((5 == parameters.m_nNonNeutrons) && (5 == parameters.m_nProtons))
            return BEAM_PARTICLE_P_P_P_P_P_P;

        return BEAM_PARTICLE_P_COMPLEX;
    }
    else if (13 == std::fabs(targetParentId))
    {
        if (0 == parameters.m_nNonNeutrons)
            return BEAM_PARTICLE_MU;
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nElectrons))
            return BEAM_PARTICLE_MU_E;
        return BEAM_PARTICLE_MU_COMPLEX;
    }
    else if (321 == targetParentId)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons))
            return BEAM_PARTICLE_KAON_PLUS_MU;
        if ((parameters.m_nNonNeutrons >= 2) && (1 == parameters.m_nKaonPlus) && (1 == parameters.m_nKaon0L))
            return BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_KAON0L_COMPLEX;
        if ((parameters.m_nNonNeutrons > 1) && (1 == parameters.m_nKaonPlus))
            return BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_COMPLEX;
        return BEAM_PARTICLE_KAON_PLUS_COMPLEX;
    }
    else if (-321 == targetParentId)
    {
        if ((1 == parameters.m_nNonNeutrons) && (1 == parameters.m_nMuons))
            return BEAM_PARTICLE_KAON_MINUS_MU;
        if ((parameters.m_nNonNeutrons >= 2) && (1 == parameters.m_nKaonMinus) && (1 == parameters.m_nKaon0L))
            return BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_KAON0L_COMPLEX;
        if ((parameters.m_nNonNeutrons > 1) && (1 == parameters.m_nKaonMinus))
            return BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_COMPLEX;
        return BEAM_PARTICLE_KAON_MINUS_COMPLEX;
    }
    else if (11 == std::fabs(targetParentId))
    {
        if (0 == parameters.m_nNonNeutrons)
            return BEAM_PARTICLE_E;
        return BEAM_PARTICLE_E_COMPLEX;
    }

    if (parameters.m_nNonNeutrons > 5)
        return BEAM_PARTICLE_COMPLEX_HIERARCHY;

    return BEAM_PARTICLE_UNKNOWN_HIERARCHY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArInteractionTypeHelper::SetInteractionParameters(const MCParticleList &mcPrimaryList, InteractionParameters &parameters)
{
    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (2112 != pMCPrimary->GetParticleId())
            ++parameters.m_nNonNeutrons;
        if (13 == std::fabs(pMCPrimary->GetParticleId()))
            ++parameters.m_nMuons;
        if (11 == std::fabs(pMCPrimary->GetParticleId()))
            ++parameters.m_nElectrons;
        else if (2212 == std::fabs(pMCPrimary->GetParticleId()))
            ++parameters.m_nProtons;
        else if (22 == pMCPrimary->GetParticleId())
            ++parameters.m_nPhotons;
        else if (211 == pMCPrimary->GetParticleId())
            ++parameters.m_nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId())
            ++parameters.m_nPiMinus;
        else if (321 == pMCPrimary->GetParticleId())
            ++parameters.m_nKaonPlus;
        else if (-321 == pMCPrimary->GetParticleId())
            ++parameters.m_nKaonMinus;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType LArInteractionTypeHelper::CosmicRayHypothesis(
    const MCParticleList &mcPrimaryList, const InteractionParameters &parameters)
{
    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsCosmicRay(mcPrimaryList.front()))
    {
        if (1 == parameters.m_nMuons)
            return COSMIC_RAY_MU;
        if (1 == parameters.m_nProtons)
            return COSMIC_RAY_P;
        if (1 == parameters.m_nElectrons)
            return COSMIC_RAY_E;
        if (1 == parameters.m_nPhotons)
            return COSMIC_RAY_PHOTON;
        else
            return COSMIC_RAY_OTHER;
    }

    return OTHER_INTERACTION;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArInteractionTypeHelper::ToString(const InteractionType interactionType)
{
    switch (interactionType)
    {
        INTERACTION_TYPE_TABLE(GET_INTERACTION_TYPE_NAME_SWITCH)
        default:
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionParameters::InteractionParameters() :
    m_nNonNeutrons(0),
    m_nMuons(0),
    m_nElectrons(0),
    m_nPhotons(0),
    m_nProtons(0),
    m_nPiPlus(0),
    m_nPiMinus(0),
    m_nPiZero(0),
    m_nKaonPlus(0),
    m_nKaonMinus(0),
    m_nKaon0L(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const int InteractionDescriptor::CC = 8192;
const int InteractionDescriptor::NC = 4096;
const int InteractionDescriptor::QE = 2560;
const int InteractionDescriptor::RES  = 2048;
const int InteractionDescriptor::DIS = 1536;
const int InteractionDescriptor::COH = 1024;
const int InteractionDescriptor::OTH = 512;
const int InteractionDescriptor::MU = 256;
const int InteractionDescriptor::E = 128;
const int InteractionDescriptor::PIZERO = 64;
const int InteractionDescriptor::PIPLUS = 32;
const int InteractionDescriptor::PIMINUS = 16;
const int InteractionDescriptor::PHOTON = 8;
const int InteractionDescriptor::NP = 6;

InteractionDescriptor::InteractionDescriptor(const bool isCC, const bool isQE, const bool isRes, const bool isDIS, const bool isCoherent,
    const bool isNumu, const bool isNue, const unsigned int nPiPlus, const unsigned int nPiMinus, const unsigned int nPhotons,
    const unsigned int nProtons) :
    m_isCC{isCC},
    m_isQE{isQE},
    m_isResonant{isRes},
    m_isDIS{isDIS},
    m_isCoherent{isCoherent},
    m_isNumu{isNumu},
    m_isNue{isNue},
    m_nPiZero{nPhotons / 2},
    m_nPiPlus{nPiPlus},
    m_nPiMinus{nPiMinus},
    m_nPhotons{nPhotons % 2},
    m_nProtons{nProtons},
    m_id{0},
    m_descriptor{}
{
    m_id += m_isCC ? CC : NC;
    m_descriptor += m_isCC ? "CC" : "NC";
    if (isQE)
    {
        m_id += QE;
        m_descriptor += "QE";
    }
    else if (isRes)
    {
        m_id += RES;
        m_descriptor += "RES";
    }
    else if (isDIS)
    {
        m_id += DIS;
        m_descriptor += "DIS";
    }
    else if (isCoherent)
    {
        m_id += COH;
        m_descriptor += "COH";
    }
    else
    {
        m_id += OTH;
        m_descriptor += "OTH";
    }

    if (m_isNumu)
    {
        m_id += MU;
        m_descriptor += "_MU";
    }
    else if (m_isNue)
    {
        m_id += E;
        m_descriptor += "_E";
    }

    m_id += m_nPiZero == 1 ? PIZERO : 0;
    m_descriptor += m_nPiZero == 1 ? "_PIZERO" : m_nPiZero > 0 ? "_NPIZERO" : "";
    m_id += m_nPiPlus == 1 ? PIPLUS : 0;
    m_descriptor += m_nPiPlus == 1 ? "_PIPLUS" : m_nPiPlus > 0 ? "_NPIPLUS" : "";
    m_id += m_nPiMinus == 1 ? PIMINUS : 0;
    m_descriptor += m_nPiMinus == 1 ? "_PIMINUS" : m_nPiMinus > 0 ? "_NPIMINUS" :  "";
    m_id += m_nPhotons == 1 ? PHOTON : 0;
    m_descriptor += m_nPhotons == 1 ? "_PHOTON" : "";
    m_id += m_nProtons < NP ? m_nProtons : NP;
    m_descriptor += m_nProtons >= NP ? "NP": m_nProtons > 0 ? "_" + std::to_string(m_nProtons) + "P" : "";
}

} // namespace lar_content
