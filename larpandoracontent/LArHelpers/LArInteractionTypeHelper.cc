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

    unsigned int nNonNeutrons(0), nMuons(0), nElectrons(0), nProtons(0), nPiPlus(0), nPiMinus(0), nPhotons(0), nKaonPlus(0), nKaonMinus(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (2112 != pMCPrimary->GetParticleId()) ++nNonNeutrons;
        if (13 == std::fabs(pMCPrimary->GetParticleId())) ++nMuons;
        if (11 == std::fabs(pMCPrimary->GetParticleId())) ++nElectrons;
        else if (2212 == std::fabs(pMCPrimary->GetParticleId())) ++nProtons;
        else if (22 == pMCPrimary->GetParticleId()) ++nPhotons;
        else if (211 == pMCPrimary->GetParticleId()) ++nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId()) ++nPiMinus;
        else if (321 == pMCPrimary->GetParticleId()) ++nKaonPlus;
        else if (-321 == pMCPrimary->GetParticleId()) ++nKaonMinus;
    }

    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsCosmicRay(mcPrimaryList.front()))
    {
        if (1 == nMuons) return COSMIC_RAY_MU;
        if (1 == nProtons) return COSMIC_RAY_P;
        if (1 == nElectrons) return COSMIC_RAY_E;
        if (1 == nPhotons) return COSMIC_RAY_PHOTON;
        else return COSMIC_RAY_OTHER;
    }

    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsBeamParticle(mcPrimaryList.front()))
    {
        if (1 == nMuons) return BEAM_PARTICLE_MU;
        if (1 == nProtons) return BEAM_PARTICLE_P;
        if (1 == nElectrons) return BEAM_PARTICLE_E;
        if (1 == nPhotons) return BEAM_PARTICLE_PHOTON;
        if (1 == nPiPlus) return BEAM_PARTICLE_PI_PLUS;
        if (1 == nPiMinus) return BEAM_PARTICLE_PI_MINUS;
        if (1 == nKaonPlus) return BEAM_PARTICLE_KAON_PLUS;
        if (1 == nKaonMinus) return BEAM_PARTICLE_KAON_MINUS;
        else return BEAM_PARTICLE_OTHER;
    }

    const MCParticle *pMCNeutrino(nullptr);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (!LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) || (pMCNeutrino && (pMCNeutrino != LArMCParticleHelper::GetParentMCParticle(pMCPrimary))))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pMCNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCPrimary);
    }

    if (!pMCNeutrino)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const int nuNuanceCode(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));

    if (1001 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCQEL_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCQEL_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCQEL_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCQEL_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCQEL_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCQEL_MU_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons)) return CCQEL_E;
        if ((2 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons)) return CCQEL_E_P;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons)) return CCQEL_E_P_P;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons)) return CCQEL_E_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons)) return CCQEL_E_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons)) return CCQEL_E_P_P_P_P_P;
    }

    if (1002 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCQEL_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCQEL_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCQEL_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCQEL_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCQEL_P_P_P_P_P;
    }

    if ((nuNuanceCode >= 1003) && (nuNuanceCode <= 1005))
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCRES_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCRES_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCRES_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCRES_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCRES_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCRES_MU_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPiPlus)) return CCRES_MU_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPhotons)) return CCRES_MU_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (2 == nPhotons)) return CCRES_MU_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_P_P_PIZERO;

        if ((1 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons)) return CCRES_E;
        if ((2 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons)) return CCRES_E_P;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons)) return CCRES_E_P_P;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons)) return CCRES_E_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons)) return CCRES_E_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons)) return CCRES_E_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (1 == nPiPlus)) return CCRES_E_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (1 == nPhotons)) return CCRES_E_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (1 == nPhotons)) return CCRES_E_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (2 == nPhotons)) return CCRES_E_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (2 == nPhotons)) return CCRES_E_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_P_P_PIZERO;
    }

    if ((nuNuanceCode >= 1006) && (nuNuanceCode <= 1009))
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCRES_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCRES_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCRES_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCRES_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCRES_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiPlus)) return NCRES_PIPLUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiPlus)) return NCRES_P_PIPLUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_PIPLUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_P_P_PIPLUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiMinus)) return NCRES_PIMINUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiMinus)) return NCRES_P_PIMINUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_PIMINUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_PIMINUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_P_PIMINUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_P_P_PIMINUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPhotons)) return NCRES_PHOTON;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPhotons)) return NCRES_P_PHOTON;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPhotons)) return NCRES_P_P_PHOTON;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_P_P_PHOTON;

        if ((2 == nNonNeutrons) && (0 == nProtons) && (2 == nPhotons)) return NCRES_PIZERO;
        if ((3 == nNonNeutrons) && (1 == nProtons) && (2 == nPhotons)) return NCRES_P_PIZERO;
        if ((4 == nNonNeutrons) && (2 == nProtons) && (2 == nPhotons)) return NCRES_P_P_PIZERO;
        if ((5 == nNonNeutrons) && (3 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (4 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (5 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_P_P_PIZERO;
    }

    if (1091 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCDIS_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCDIS_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCDIS_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCDIS_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCDIS_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCDIS_MU_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPiPlus)) return CCDIS_MU_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPhotons)) return CCDIS_MU_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPhotons)) return CCDIS_MU_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (2 == nPhotons)) return CCDIS_MU_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (2 == nPhotons)) return CCDIS_MU_P_P_P_P_P_PIZERO;
    }

    if (1092 == nuNuanceCode)
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCDIS_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCDIS_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCDIS_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCDIS_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCDIS_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiPlus)) return NCDIS_PIPLUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiPlus)) return NCDIS_P_PIPLUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_PIPLUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiPlus)) return NCDIS_P_P_P_P_P_PIPLUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiMinus)) return NCDIS_PIMINUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiMinus)) return NCDIS_P_PIMINUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_PIMINUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_P_PIMINUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_P_P_PIMINUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiMinus)) return NCDIS_P_P_P_P_P_PIMINUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPhotons)) return NCDIS_PHOTON;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPhotons)) return NCDIS_P_PHOTON;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_PHOTON;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPhotons)) return NCDIS_P_P_P_P_P_PHOTON;

        if ((2 == nNonNeutrons) && (0 == nProtons) && (2 == nPhotons)) return NCDIS_PIZERO;
        if ((3 == nNonNeutrons) && (1 == nProtons) && (2 == nPhotons)) return NCDIS_P_PIZERO;
        if ((4 == nNonNeutrons) && (2 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_PIZERO;
        if ((5 == nNonNeutrons) && (3 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (4 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (5 == nProtons) && (2 == nPhotons)) return NCDIS_P_P_P_P_P_PIZERO;
    }

    if (1096 == nuNuanceCode) return NCCOH;
    if (1097 == nuNuanceCode) return CCCOH;

    return OTHER_INTERACTION;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArInteractionTypeHelper::InteractionType LArInteractionTypeHelper::GetTestBeamHierarchyInteractionType(const MCParticleList &mcPrimaryList)
{
    if (mcPrimaryList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    unsigned int nNonNeutrons(0), nMuons(0), nElectrons(0), nProtons(0), nPiPlus(0), nPiMinus(0), nPhotons(0), nKaonPlus(0), nKaonMinus(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (2112 != pMCPrimary->GetParticleId()) ++nNonNeutrons;
        if (13 == std::fabs(pMCPrimary->GetParticleId())) ++nMuons;
        if (11 == std::fabs(pMCPrimary->GetParticleId())) ++nElectrons;
        else if (2212 == std::fabs(pMCPrimary->GetParticleId())) ++nProtons;
        else if (22 == pMCPrimary->GetParticleId()) ++nPhotons;
        else if (211 == pMCPrimary->GetParticleId()) ++nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId()) ++nPiMinus;
        else if (321 == pMCPrimary->GetParticleId()) ++nKaonPlus;
        else if (-321 == pMCPrimary->GetParticleId()) ++nKaonMinus;
    }

    if ((1 == mcPrimaryList.size()) && LArMCParticleHelper::IsCosmicRay(mcPrimaryList.front()))
    {
        if (1 == nMuons) return COSMIC_RAY_MU;
        if (1 == nProtons) return COSMIC_RAY_P;
        if (1 == nElectrons) return COSMIC_RAY_E;
        if (1 == nPhotons) return COSMIC_RAY_PHOTON;
        else return COSMIC_RAY_OTHER;
    }

    int targetParentId(0);

    for (const MCParticle *const pMCPrimary : mcPrimaryList)
    {
        if (LArMCParticleHelper::GetParentMCParticle(pMCPrimary) != pMCPrimary)
            continue;

        if (13 == std::fabs(pMCPrimary->GetParticleId())) --nMuons;
        else if (11 == std::fabs(pMCPrimary->GetParticleId())) --nElectrons;
        else if (2212 == std::fabs(pMCPrimary->GetParticleId())) --nProtons;
        else if (22 == pMCPrimary->GetParticleId()) --nPhotons;
        else if (211 == pMCPrimary->GetParticleId()) --nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId()) --nPiMinus;
        else if (321 == pMCPrimary->GetParticleId()) --nKaonPlus;
        else if (-321 == pMCPrimary->GetParticleId()) --nKaonMinus;

        --nNonNeutrons;
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

    const int nPiZero(piZero.size());
    const int nKaon0L(kaon0L.size());

    if (211 == targetParentId)
    {
        if ((1 == nNonNeutrons) && (1 == nPiPlus)) return BEAM_PARTICLE_PI_PLUS_PI_PLUS;
        if ((2 == nNonNeutrons) && (1 == nPiPlus) && (1 == nPhotons)) return BEAM_PARTICLE_PI_PLUS_PI_PLUS_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nPiPlus) && (1 == nPiZero)) return BEAM_PARTICLE_PI_PLUS_PI_PLUS_PIZERO;
        return BEAM_PARTICLE_PI_PLUS_COMPLEX;
    }
    else if (-211 == targetParentId)
    {
        if ((1 == nNonNeutrons) && (1 == nPiMinus)) return BEAM_PARTICLE_PI_MINUS_PI_MINUS;
        if ((2 == nNonNeutrons) && (1 == nPiMinus) && (1 == nPhotons)) return BEAM_PARTICLE_PI_MINUS_PI_MINUS_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nPiMinus) && (1 == nPiZero)) return BEAM_PARTICLE_PI_MINUS_PI_MINUS_PIZERO;
        return BEAM_PARTICLE_PI_MINUS_COMPLEX;
    }
    else if (2212 == targetParentId)
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return BEAM_PARTICLE_P_P;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPhotons)) return BEAM_PARTICLE_P_P_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nProtons) && (2 == nPhotons)) return BEAM_PARTICLE_P_P_PHOTON_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nProtons) && (3 == nPhotons)) return BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nProtons) && (4 == nPhotons)) return BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON_PHOTON;

        if ((2 == nNonNeutrons) && (2 == nProtons)) return BEAM_PARTICLE_P_P_P;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPhotons)) return BEAM_PARTICLE_P_P_P_PHOTON;
        if ((4 == nNonNeutrons) && (2 == nProtons) && (2 == nPhotons)) return BEAM_PARTICLE_P_P_P_PHOTON_PHOTON;
        if ((5 == nNonNeutrons) && (2 == nProtons) && (3 == nPhotons)) return BEAM_PARTICLE_P_P_P_PHOTON_PHOTON_PHOTON;

        if ((3 == nNonNeutrons) && (3 == nProtons)) return BEAM_PARTICLE_P_P_P_P;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPhotons)) return BEAM_PARTICLE_P_P_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (3 == nProtons) && (2 == nPhotons)) return BEAM_PARTICLE_P_P_P_P_PHOTON_PHOTON;

        if ((4 == nNonNeutrons) && (4 == nProtons)) return BEAM_PARTICLE_P_P_P_P_P;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPhotons)) return BEAM_PARTICLE_P_P_P_P_P_PHOTON;

        if ((5 == nNonNeutrons) && (5 == nProtons)) return BEAM_PARTICLE_P_P_P_P_P_P;

        return BEAM_PARTICLE_P_COMPLEX;
    }
    else if (13 == std::fabs(targetParentId))
    {
        if (0 == nNonNeutrons) return BEAM_PARTICLE_MU;
        if ((1 == nNonNeutrons) && (1 == nElectrons)) return BEAM_PARTICLE_MU_E;
        return BEAM_PARTICLE_MU_COMPLEX;
    }
    else if (321 == targetParentId)
    {
        if ((1 == nNonNeutrons) && (1 == nMuons)) return BEAM_PARTICLE_KAON_PLUS_MU;
        if ((nNonNeutrons >= 2) && (1 == nKaonPlus) && (1 == nKaon0L)) return BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_KAON0L_COMPLEX;
        if ((nNonNeutrons > 1) && (1 == nKaonPlus)) return BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_COMPLEX;
        return BEAM_PARTICLE_KAON_PLUS_COMPLEX;
    }
    else if (-321 == targetParentId)
    {
        if ((1 == nNonNeutrons) && (1 == nMuons)) return BEAM_PARTICLE_KAON_MINUS_MU;
        if ((nNonNeutrons >= 2) && (1 == nKaonMinus) && (1 == nKaon0L)) return BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_KAON0L_COMPLEX;
        if ((nNonNeutrons > 1) && (1 == nKaonMinus)) return BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_COMPLEX;
        return BEAM_PARTICLE_KAON_MINUS_COMPLEX;
    }
    else if (11 == std::fabs(targetParentId))
    {
        if (0 == nNonNeutrons) return BEAM_PARTICLE_E;
        return BEAM_PARTICLE_E_COMPLEX;
    }

    if (nNonNeutrons > 5) return BEAM_PARTICLE_COMPLEX_HIERARCHY;

    return BEAM_PARTICLE_UNKNOWN_HIERARCHY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArInteractionTypeHelper::ToString(const InteractionType interactionType)
{
    switch (interactionType)
    {
    case CCQEL_MU: return "CCQEL_MU";
    case CCQEL_MU_P: return "CCQEL_MU_P";
    case CCQEL_MU_P_P: return "CCQEL_MU_P_P";
    case CCQEL_MU_P_P_P: return "CCQEL_MU_P_P_P";
    case CCQEL_MU_P_P_P_P: return "CCQEL_MU_P_P_P_P";
    case CCQEL_MU_P_P_P_P_P: return "CCQEL_MU_P_P_P_P_P";
    case CCQEL_E: return "CCQEL_E";
    case CCQEL_E_P: return "CCQEL_E_P";
    case CCQEL_E_P_P: return "CCQEL_E_P_P";
    case CCQEL_E_P_P_P: return "CCQEL_E_P_P_P";
    case CCQEL_E_P_P_P_P: return "CCQEL_E_P_P_P_P";
    case CCQEL_E_P_P_P_P_P: return "CCQEL_E_P_P_P_P_P";
    case NCQEL_P: return "NCQEL_P";
    case NCQEL_P_P: return "NCQEL_P_P";
    case NCQEL_P_P_P: return "NCQEL_P_P_P";
    case NCQEL_P_P_P_P: return "NCQEL_P_P_P_P";
    case NCQEL_P_P_P_P_P: return "NCQEL_P_P_P_P_P";
    case CCRES_MU: return "CCRES_MU";
    case CCRES_MU_P: return "CCRES_MU_P";
    case CCRES_MU_P_P: return "CCRES_MU_P_P";
    case CCRES_MU_P_P_P: return "CCRES_MU_P_P_P";
    case CCRES_MU_P_P_P_P: return "CCRES_MU_P_P_P_P";
    case CCRES_MU_P_P_P_P_P: return "CCRES_MU_P_P_P_P_P";
    case CCRES_MU_PIPLUS: return "CCRES_MU_PIPLUS";
    case CCRES_MU_P_PIPLUS: return "CCRES_MU_P_PIPLUS";
    case CCRES_MU_P_P_PIPLUS: return "CCRES_MU_P_P_PIPLUS";
    case CCRES_MU_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_P_PIPLUS";
    case CCRES_MU_PHOTON: return "CCRES_MU_PHOTON";
    case CCRES_MU_P_PHOTON: return "CCRES_MU_P_PHOTON";
    case CCRES_MU_P_P_PHOTON: return "CCRES_MU_P_P_PHOTON";
    case CCRES_MU_P_P_P_PHOTON: return "CCRES_MU_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_P_PHOTON";
    case CCRES_MU_PIZERO: return "CCRES_MU_PIZERO";
    case CCRES_MU_P_PIZERO: return "CCRES_MU_P_PIZERO";
    case CCRES_MU_P_P_PIZERO: return "CCRES_MU_P_P_PIZERO";
    case CCRES_MU_P_P_P_PIZERO: return "CCRES_MU_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_P_PIZERO";
    case CCRES_E: return "CCRES_E";
    case CCRES_E_P: return "CCRES_E_P";
    case CCRES_E_P_P: return "CCRES_E_P_P";
    case CCRES_E_P_P_P: return "CCRES_E_P_P_P";
    case CCRES_E_P_P_P_P: return "CCRES_E_P_P_P_P";
    case CCRES_E_P_P_P_P_P: return "CCRES_E_P_P_P_P_P";
    case CCRES_E_PIPLUS: return "CCRES_E_PIPLUS";
    case CCRES_E_P_PIPLUS: return "CCRES_E_P_PIPLUS";
    case CCRES_E_P_P_PIPLUS: return "CCRES_E_P_P_PIPLUS";
    case CCRES_E_P_P_P_PIPLUS: return "CCRES_E_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_P_PIPLUS";
    case CCRES_E_PHOTON: return "CCRES_E_PHOTON";
    case CCRES_E_P_PHOTON: return "CCRES_E_P_PHOTON";
    case CCRES_E_P_P_PHOTON: return "CCRES_E_P_P_PHOTON";
    case CCRES_E_P_P_P_PHOTON: return "CCRES_E_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_P_PHOTON";
    case CCRES_E_PIZERO: return "CCRES_E_PIZERO";
    case CCRES_E_P_PIZERO: return "CCRES_E_P_PIZERO";
    case CCRES_E_P_P_PIZERO: return "CCRES_E_P_P_PIZERO";
    case CCRES_E_P_P_P_PIZERO: return "CCRES_E_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_P_PIZERO";
    case NCRES_P: return "NCRES_P";
    case NCRES_P_P: return "NCRES_P_P";
    case NCRES_P_P_P: return "NCRES_P_P_P";
    case NCRES_P_P_P_P: return "NCRES_P_P_P_P";
    case NCRES_P_P_P_P_P: return "NCRES_P_P_P_P_P";
    case NCRES_PIPLUS: return "NCRES_PIPLUS";
    case NCRES_P_PIPLUS: return "NCRES_P_PIPLUS";
    case NCRES_P_P_PIPLUS: return "NCRES_P_P_PIPLUS";
    case NCRES_P_P_P_PIPLUS: return "NCRES_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_P_PIPLUS";
    case NCRES_PIMINUS: return "NCRES_PIMINUS";
    case NCRES_P_PIMINUS: return "NCRES_P_PIMINUS";
    case NCRES_P_P_PIMINUS: return "NCRES_P_P_PIMINUS";
    case NCRES_P_P_P_PIMINUS: return "NCRES_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_P_PIMINUS";
    case NCRES_PHOTON: return "NCRES_PHOTON";
    case NCRES_P_PHOTON: return "NCRES_P_PHOTON";
    case NCRES_P_P_PHOTON: return "NCRES_P_P_PHOTON";
    case NCRES_P_P_P_PHOTON: return "NCRES_P_P_P_PHOTON";
    case NCRES_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_PHOTON";
    case NCRES_P_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_P_PHOTON";
    case NCRES_PIZERO: return "NCRES_PIZERO";
    case NCRES_P_PIZERO: return "NCRES_P_PIZERO";
    case NCRES_P_P_PIZERO: return "NCRES_P_P_PIZERO";
    case NCRES_P_P_P_PIZERO: return "NCRES_P_P_P_PIZERO";
    case NCRES_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_PIZERO";
    case NCRES_P_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_P_PIZERO";
    case CCDIS_MU: return "CCDIS_MU";
    case CCDIS_MU_P: return "CCDIS_MU_P";
    case CCDIS_MU_P_P: return "CCDIS_MU_P_P";
    case CCDIS_MU_P_P_P: return "CCDIS_MU_P_P_P";
    case CCDIS_MU_P_P_P_P: return "CCDIS_MU_P_P_P_P";
    case CCDIS_MU_P_P_P_P_P: return "CCDIS_MU_P_P_P_P_P";
    case CCDIS_MU_PIPLUS: return "CCDIS_MU_PIPLUS";
    case CCDIS_MU_P_PIPLUS: return "CCDIS_MU_P_PIPLUS";
    case CCDIS_MU_P_P_PIPLUS: return "CCDIS_MU_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_P_PIPLUS";
    case CCDIS_MU_P_P_P_P_P_PIPLUS: return "CCDIS_MU_P_P_P_P_P_PIPLUS";
    case CCDIS_MU_PHOTON: return "CCDIS_MU_PHOTON";
    case CCDIS_MU_P_PHOTON: return "CCDIS_MU_P_PHOTON";
    case CCDIS_MU_P_P_PHOTON: return "CCDIS_MU_P_P_PHOTON";
    case CCDIS_MU_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_PHOTON";
    case CCDIS_MU_P_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_P_PHOTON";
    case CCDIS_MU_P_P_P_P_P_PHOTON: return "CCDIS_MU_P_P_P_P_P_PHOTON";
    case CCDIS_MU_PIZERO: return "CCDIS_MU_PIZERO";
    case CCDIS_MU_P_PIZERO: return "CCDIS_MU_P_PIZERO";
    case CCDIS_MU_P_P_PIZERO: return "CCDIS_MU_P_P_PIZERO";
    case CCDIS_MU_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_PIZERO";
    case CCDIS_MU_P_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_P_PIZERO";
    case CCDIS_MU_P_P_P_P_P_PIZERO: return "CCDIS_MU_P_P_P_P_P_PIZERO";
    case NCDIS_P: return "NCDIS_P";
    case NCDIS_P_P: return "NCDIS_P_P";
    case NCDIS_P_P_P: return "NCDIS_P_P_P";
    case NCDIS_P_P_P_P: return "NCDIS_P_P_P_P";
    case NCDIS_P_P_P_P_P: return "NCDIS_P_P_P_P_P";
    case NCDIS_PIPLUS: return "NCDIS_PIPLUS";
    case NCDIS_P_PIPLUS: return "NCDIS_P_PIPLUS";
    case NCDIS_P_P_PIPLUS: return "NCDIS_P_P_PIPLUS";
    case NCDIS_P_P_P_PIPLUS: return "NCDIS_P_P_P_PIPLUS";
    case NCDIS_P_P_P_P_PIPLUS: return "NCDIS_P_P_P_P_PIPLUS";
    case NCDIS_P_P_P_P_P_PIPLUS: return "NCDIS_P_P_P_P_P_PIPLUS";
    case NCDIS_PIMINUS: return "NCDIS_PIMINUS";
    case NCDIS_P_PIMINUS: return "NCDIS_P_PIMINUS";
    case NCDIS_P_P_PIMINUS: return "NCDIS_P_P_PIMINUS";
    case NCDIS_P_P_P_PIMINUS: return "NCDIS_P_P_P_PIMINUS";
    case NCDIS_P_P_P_P_PIMINUS: return "NCDIS_P_P_P_P_PIMINUS";
    case NCDIS_P_P_P_P_P_PIMINUS: return "NCDIS_P_P_P_P_P_PIMINUS";
    case NCDIS_PHOTON: return "NCDIS_PHOTON";
    case NCDIS_P_PHOTON: return "NCDIS_P_PHOTON";
    case NCDIS_P_P_PHOTON: return "NCDIS_P_P_PHOTON";
    case NCDIS_P_P_P_PHOTON: return "NCDIS_P_P_P_PHOTON";
    case NCDIS_P_P_P_P_PHOTON: return "NCDIS_P_P_P_P_PHOTON";
    case NCDIS_P_P_P_P_P_PHOTON: return "NCDIS_P_P_P_P_P_PHOTON";
    case NCDIS_PIZERO: return "NCDIS_PIZERO";
    case NCDIS_P_PIZERO: return "NCDIS_P_PIZERO";
    case NCDIS_P_P_PIZERO: return "NCDIS_P_P_PIZERO";
    case NCDIS_P_P_P_PIZERO: return "NCDIS_P_P_P_PIZERO";
    case NCDIS_P_P_P_P_PIZERO: return "NCDIS_P_P_P_P_PIZERO";
    case NCDIS_P_P_P_P_P_PIZERO: return "NCDIS_P_P_P_P_P_PIZERO";
    case CCCOH: return "CCCOH";
    case NCCOH: return "NCCOH";
    case COSMIC_RAY_MU: return "COSMIC_RAY_MU";
    case COSMIC_RAY_P: return "COSMIC_RAY_P";
    case COSMIC_RAY_E: return "COSMIC_RAY_E";
    case COSMIC_RAY_PHOTON: return "COSMIC_RAY_PHOTON";
    case COSMIC_RAY_OTHER: return "COSMIC_RAY_OTHER";
    case BEAM_PARTICLE_MU: return "BEAM_PARTICLE_MU";
    case BEAM_PARTICLE_P: return "BEAM_PARTICLE_P";
    case BEAM_PARTICLE_E: return "BEAM_PARTICLE_E";
    case BEAM_PARTICLE_PHOTON: return "BEAM_PARTICLE_PHOTON";
    case BEAM_PARTICLE_PI_PLUS: return "BEAM_PARTICLE_PI_PLUS";
    case BEAM_PARTICLE_PI_MINUS: return "BEAM_PARTICLE_PI_MINUS";
    case BEAM_PARTICLE_KAON_PLUS: return "BEAM_PARTICLE_KAON_PLUS";
    case BEAM_PARTICLE_KAON_MINUS: return "BEAM_PARTICLE_KAON_MINUS";
    case BEAM_PARTICLE_OTHER: return "BEAM_PARTICLE_OTHER";
    case BEAM_PARTICLE_PI_PLUS_PI_PLUS: return "BEAM_PARTICLE_PI_PLUS_PI_PLUS";
    case BEAM_PARTICLE_PI_PLUS_PI_PLUS_PHOTON: return "BEAM_PARTICLE_PI_PLUS_PI_PLUS_PHOTON";
    case BEAM_PARTICLE_PI_PLUS_PI_PLUS_PIZERO: return "BEAM_PARTICLE_PI_PLUS_PI_PLUS_PIZERO";
    case BEAM_PARTICLE_PI_PLUS_COMPLEX: return "BEAM_PARTICLE_PI_PLUS_COMPLEX";
    case BEAM_PARTICLE_PI_MINUS_PI_MINUS: return "BEAM_PARTICLE_PI_MINUS_PI_MINUS";
    case BEAM_PARTICLE_PI_MINUS_PI_MINUS_PHOTON: return "BEAM_PARTICLE_PI_MINUS_PI_MINUS_PHOTON";
    case BEAM_PARTICLE_PI_MINUS_PI_MINUS_PIZERO: return "BEAM_PARTICLE_PI_MINUS_PI_MINUS_PIZERO";
    case BEAM_PARTICLE_PI_MINUS_COMPLEX: return "BEAM_PARTICLE_PI_MINUS_COMPLEX";
    case BEAM_PARTICLE_P_P: return "BEAM_PARTICLE_P_P";
    case BEAM_PARTICLE_P_P_PHOTON: return "BEAM_PARTICLE_P_P_PHOTON";
    case BEAM_PARTICLE_P_P_PHOTON_PHOTON: return "BEAM_PARTICLE_P_P_PHOTON_PHOTON";
    case BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON: return "BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON";
    case BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON_PHOTON: return "BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON_PHOTON";
    case BEAM_PARTICLE_P_P_P: return "BEAM_PARTICLE_P_P_P";
    case BEAM_PARTICLE_P_P_P_PHOTON: return "BEAM_PARTICLE_P_P_P_PHOTON";
    case BEAM_PARTICLE_P_P_P_PHOTON_PHOTON: return "BEAM_PARTICLE_P_P_P_PHOTON_PHOTON";
    case BEAM_PARTICLE_P_P_P_PHOTON_PHOTON_PHOTON: return "BEAM_PARTICLE_P_P_P_PHOTON_PHOTON_PHOTON";
    case BEAM_PARTICLE_P_P_P_P: return "BEAM_PARTICLE_P_P_P_P";
    case BEAM_PARTICLE_P_P_P_P_PHOTON: return "BEAM_PARTICLE_P_P_P_P_PHOTON";
    case BEAM_PARTICLE_P_P_P_P_PHOTON_PHOTON: return "BEAM_PARTICLE_P_P_P_P_PHOTON_PHOTON";
    case BEAM_PARTICLE_P_P_P_P_P: return "BEAM_PARTICLE_P_P_P_P_P";
    case BEAM_PARTICLE_P_P_P_P_P_PHOTON: return "BEAM_PARTICLE_P_P_P_P_P_PHOTON";
    case BEAM_PARTICLE_P_P_P_P_P_P: return "BEAM_PARTICLE_P_P_P_P_P_P";
    case BEAM_PARTICLE_P_COMPLEX: return "BEAM_PARTICLE_P_COMPLEX";
    case BEAM_PARTICLE_MU_E: return "BEAM_PARTICLE_MU_E";
    case BEAM_PARTICLE_MU_COMPLEX: return "BEAM_PARTICLE_MU_COMPLEX";
    case BEAM_PARTICLE_KAON_PLUS_MU: return "BEAM_PARTICLE_KAON_PLUS_MU";
    case BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_KAON0L_COMPLEX: return "BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_KAON0L_COMPLEX";
    case BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_COMPLEX: return "BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_COMPLEX";
    case BEAM_PARTICLE_KAON_PLUS_COMPLEX: return "BEAM_PARTICLE_KAON_PLUS_COMPLEX";
    case BEAM_PARTICLE_KAON_MINUS_MU: return "BEAM_PARTICLE_KAON_MINUS_MU";
    case BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_KAON0L_COMPLEX: return "BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_KAON0L_COMPLEX";
    case BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_COMPLEX: return "BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_COMPLEX";
    case BEAM_PARTICLE_KAON_MINUS_COMPLEX: return "BEAM_PARTICLE_KAON_MINUS_COMPLEX";
    case BEAM_PARTICLE_E_COMPLEX: return "BEAM_PARTICLE_E_COMPLEX";
    case BEAM_PARTICLE_COMPLEX_HIERARCHY: return "BEAM_PARTICLE_COMPLEX_HIERARCHY";
    case BEAM_PARTICLE_UNKNOWN_HIERARCHY: return "BEAM_PARTICLE_UNKNOWN_HIERARCHY";

    case OTHER_INTERACTION: return "OTHER_INTERACTION";
    case ALL_INTERACTIONS: return "ALL_INTERACTIONS";
    default: return "UNKNOWN";
    }
}

} // namespace lar_content
