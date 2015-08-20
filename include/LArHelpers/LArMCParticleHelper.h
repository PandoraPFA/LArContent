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

#include <unordered_map>

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
     *  @brief  Whether a mc particle is visible (i.e. long-lived charged particle)
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return boolean
     */
    static bool IsVisible(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get neutrino MC particles from an input MC particle list
     * 
     *  @param  pMCParticleList the input MC particle list
     *  @param  trueNeutrinos to receive the list of neutrino MC particles
     */
    static void GetTrueNeutrinos(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleList &trueNeutrinos);

    /**
     *  @brief  Get the primary parent mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetPrimaryMCParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the parent mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the parent mc particle
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

    typedef std::unordered_map<const pandora::MCParticle*, const pandora::MCParticle*> MCRelationMap;

    /**
     *  @brief  Whether a provided mc particle matches the implemented definition of being primary
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return boolean
     */
    static bool IsPrimary(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) and their primary parent mc particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcPrimaryMap the output mapping between mc particles and their parents
     */
    static void GetMCPrimaryMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap);

    /**
     *  @brief  Get vector of primary MC particles from an input list of MC particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcPrimaryVector the output mc particle vector
     */
    static void GetPrimaryMCParticleList(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &mcPrimaryVector);

    /**
     *  @brief  Get vector of primary neutrinos from an input list of MC particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcNeutrinoVector the output mc particle vector
     */
    static void GetNeutrinoMCParticleList(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &mcNeutrinoVector);

    /**
     *  @brief  Find the mc particle making the largest contribution to 2D clusters in a specified pfo
     *
     *  @param  pT address of the pfo to examine
     * 
     *  @return address of the main mc particle
     */
    static const pandora::MCParticle *GetMainMCParticle(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Find the primary mc particle making the largest contribution to 2D clusters in a specified pfo
     *
     *  @param  pT address of the pfo to examine
     *  @param  mcPrimaryMap the provided mapping between mc particles and their parents
     * 
     *  @return address of the main mc primary
     */
    static const pandora::MCParticle *GetMainMCPrimary(const pandora::ParticleFlowObject *const pPfo, const MCRelationMap &mcPrimaryMap);

    /**
     *  @brief  Sort mc particles by their source
     *
     *  @param  pLhs address of first mc particle
     *  @param  pRhs address of second mc particle
     */
    static bool SortBySource(const pandora::MCParticle *const pLhs, const pandora::MCParticle *const pRhs);

    /**
     *  @brief  Sort mc particles by their momentum
     *
     *  @param  pLhs address of first mc particle
     *  @param  pRhs address of second mc particle
     */
    static bool SortByMomentum(const pandora::MCParticle *const pLhs, const pandora::MCParticle *const pRhs);
};

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_HELPER_H
