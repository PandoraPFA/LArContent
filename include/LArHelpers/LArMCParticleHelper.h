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
     *  @brief  Whether a mc particle is a final-state particle from a neutrino or antineutrino interaction
     * 
     *  @param  pMCParticle the input mc particle
     * 
     *  @return boolean
     */
    static bool IsNeutrinoFinalState(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is produced by the interation of a neutrino or antineutrino
     * 
     *  @param  pMCParticle the input mc particle
     * 
     *  @return boolean
     */
    static bool IsNeutrinoInduced(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is a neutrino or antineutrino
     * 
     *  @param  pMCParticle the input mc particle
     * 
     *  @return boolean
     */
    static bool IsNeutrino(const pandora::MCParticle *const pMCParticle);    

    /**
     *  @brief  Get the primary parent mc particle
     * 
     *  @param  pMCParticle the input mc particle
     * 
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetParentMCParticle(const pandora::MCParticle *const pMCParticle);

     /**
     *  @brief  Get parent neutrino or antineutrino
     * 
     *  @param  pMCParticle the input mc particle
     * 
     *  @return address of primary neutrino mc particle
     */
    static const pandora::MCParticle *GetParentNeutrino(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get parent neutrino or antineutrino pdg code
     * 
     *  @param  pMCParticle the input mc particle
     * 
     *  @return pdg code of neutrino (or zero, otherwise)
     */
    static int GetParentNeutrinoId(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a cluster is the product of a neutrino interaction
     * 
     *  @param  pCluster the input cluster
     *  @param  minWeight minimum threshold for the combined weight of neutrino-induced particles associated with the cluster
     * 
     *  @return boolean
     */
    static bool IsNeutrinoInduced(const pandora::Cluster *const pCluster, const float minWeight = 0.f);

    /**
     *  @brief  Whether a hit is the product of a neutrino interaction
     * 
     *  @param  pCaloHit the input calo hit
     *  @param  minWeight minimum threshold for the combined weight of neutrino-induced particles associated with the hit
     * 
     *  @return boolean
     */
    static bool IsNeutrinoInduced(const pandora::CaloHit *const pCaloHit, const float minWeight = 0.f);

    /**
     *  @brief  Calculate the combined weight of neutrino-induced particles associated with a cluster
     * 
     *  @param  pCluster the input cluster
     * 
     *  @return the combined weight of neutrino-induced particles associated with the cluster
     */
    static float GetNeutrinoWeight(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Calculate the combined weight of neutrino-induced particles associated with a calo hit
     * 
     *  @param  pCaloHit the input calo hit
     * 
     *  @return the combined weight of neutrino-induced particles associated with the calo hit
     */
    static float GetNeutrinoWeight(const pandora::CaloHit *const pCaloHit);
};

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_HELPER_H
