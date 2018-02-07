/**
 *  @file   larpandoracontent/LArHelpers/LArInteractionTypeHelper.h
 *
 *  @brief  Header file for the interaction type helper class.
 *
 *  $Log: $
 */
#ifndef LAR_INTERACTION_TYPE_HELPER_H
#define LAR_INTERACTION_TYPE_HELPER_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"

namespace lar_content
{

/**
 *  @brief  LArInteractionTypeHelper class
 */
class LArInteractionTypeHelper
{
public:
    /**
     *  @brief   InteractionType enum
     */
    enum InteractionType : unsigned int
    {
        CCQEL_MU,
        CCQEL_MU_P,
        CCQEL_MU_P_P,
        CCQEL_MU_P_P_P,
        CCQEL_MU_P_P_P_P,
        CCQEL_MU_P_P_P_P_P,
        CCQEL_E,
        CCQEL_E_P,
        CCQEL_E_P_P,
        CCQEL_E_P_P_P,
        CCQEL_E_P_P_P_P,
        CCQEL_E_P_P_P_P_P,
        NCQEL_P,
        NCQEL_P_P,
        NCQEL_P_P_P,
        NCQEL_P_P_P_P,
        NCQEL_P_P_P_P_P,
        CCRES_MU,
        CCRES_MU_P,
        CCRES_MU_P_P,
        CCRES_MU_P_P_P,
        CCRES_MU_P_P_P_P,
        CCRES_MU_P_P_P_P_P,
        CCRES_MU_PIPLUS,
        CCRES_MU_P_PIPLUS,
        CCRES_MU_P_P_PIPLUS,
        CCRES_MU_P_P_P_PIPLUS,
        CCRES_MU_P_P_P_P_PIPLUS,
        CCRES_MU_P_P_P_P_P_PIPLUS,
        CCRES_MU_PHOTON,
        CCRES_MU_P_PHOTON,
        CCRES_MU_P_P_PHOTON,
        CCRES_MU_P_P_P_PHOTON,
        CCRES_MU_P_P_P_P_PHOTON,
        CCRES_MU_P_P_P_P_P_PHOTON,
        CCRES_MU_PIZERO,
        CCRES_MU_P_PIZERO,
        CCRES_MU_P_P_PIZERO,
        CCRES_MU_P_P_P_PIZERO,
        CCRES_MU_P_P_P_P_PIZERO,
        CCRES_MU_P_P_P_P_P_PIZERO,
        CCRES_E,
        CCRES_E_P,
        CCRES_E_P_P,
        CCRES_E_P_P_P,
        CCRES_E_P_P_P_P,
        CCRES_E_P_P_P_P_P,
        CCRES_E_PIPLUS,
        CCRES_E_P_PIPLUS,
        CCRES_E_P_P_PIPLUS,
        CCRES_E_P_P_P_PIPLUS,
        CCRES_E_P_P_P_P_PIPLUS,
        CCRES_E_P_P_P_P_P_PIPLUS,
        CCRES_E_PHOTON,
        CCRES_E_P_PHOTON,
        CCRES_E_P_P_PHOTON,
        CCRES_E_P_P_P_PHOTON,
        CCRES_E_P_P_P_P_PHOTON,
        CCRES_E_P_P_P_P_P_PHOTON,
        CCRES_E_PIZERO,
        CCRES_E_P_PIZERO,
        CCRES_E_P_P_PIZERO,
        CCRES_E_P_P_P_PIZERO,
        CCRES_E_P_P_P_P_PIZERO,
        CCRES_E_P_P_P_P_P_PIZERO,
        NCRES_P,
        NCRES_P_P,
        NCRES_P_P_P,
        NCRES_P_P_P_P,
        NCRES_P_P_P_P_P,
        NCRES_PIPLUS,
        NCRES_P_PIPLUS,
        NCRES_P_P_PIPLUS,
        NCRES_P_P_P_PIPLUS,
        NCRES_P_P_P_P_PIPLUS,
        NCRES_P_P_P_P_P_PIPLUS,
        NCRES_PIMINUS,
        NCRES_P_PIMINUS,
        NCRES_P_P_PIMINUS,
        NCRES_P_P_P_PIMINUS,
        NCRES_P_P_P_P_PIMINUS,
        NCRES_P_P_P_P_P_PIMINUS,
        NCRES_PHOTON,
        NCRES_P_PHOTON,
        NCRES_P_P_PHOTON,
        NCRES_P_P_P_PHOTON,
        NCRES_P_P_P_P_PHOTON,
        NCRES_P_P_P_P_P_PHOTON,
        NCRES_PIZERO,
        NCRES_P_PIZERO,
        NCRES_P_P_PIZERO,
        NCRES_P_P_P_PIZERO,
        NCRES_P_P_P_P_PIZERO,
        NCRES_P_P_P_P_P_PIZERO,
        CCDIS,
        NCDIS,
        CCCOH,
        NCCOH,
        COSMIC_RAY,
        BEAM_PARTICLE,
        OTHER_INTERACTION,
        ALL_INTERACTIONS
    };

    /**
     *  @brief  Get the interaction type of an event
     *
     *  @param  nuanceCode the nuance code
     *  @param  mcPrimaryList the list of primary mc particles
     *
     *  @return interaction type
     */
    static InteractionType GetInteractionType(const pandora::MCParticleList &mcPrimaryList);

    /**
     *  @brief  Get a string representation of an interaction type
     *
     *  @param  interactionType the interaction type
     *
     *  @return string
     */
    static std::string ToString(const InteractionType interactionType);
};

} // namespace lar_content

#endif // #ifndef LAR_INTERACTION_TYPE_HELPER_H

