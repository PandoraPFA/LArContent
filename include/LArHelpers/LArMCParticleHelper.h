/**
 *  @file   LArContent/include/LArHelpers/LArMCParticleHelper.h
 * 
 *  @brief  Header file for the lar monte carlo particle helper helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_MC_PARTICLE_HELPER_H
#define LAR_MC_PARTICLE_HELPER_H 1

#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
 *  @brief  LArMCParticleHelper class
 */
class LArMCParticleHelper
{
public:
    /**
     *  @brief  Whether a mc particle is a final-state particle from a neutrino (or antineutrino) interaction
     * 
     *  @param  pMCParticle
     * 
     *  @return boolean
     */
    static bool IsNeutrinoFinalState(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is a neutrino or (antineutrino)
     * 
     *  @param  pMCParticle
     * 
     *  @return boolean
     */
    static bool IsNeutrino(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the primary parent mc particle
     * 
     *  @param  pMCParticle
     * 
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetParentMCParticle(const pandora::MCParticle *const pMCParticle);

     /**
     *  @brief  Get parent neutrino or antineutrino
     * 
     *  @param  pMCParticle
     * 
     *  @return address of primary neutrino mc particle
     */
    static const pandora::MCParticle *GetParentNeutrino(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get parent neutrino or antineutrino pdg code
     * 
     *  @param  pMCParticle
     * 
     *  @return pdg code of neutrino (or zero, otherwise)
     */
    static int GetParentNeutrinoId(const pandora::MCParticle *const pMCParticle);
};

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_HELPER_H
