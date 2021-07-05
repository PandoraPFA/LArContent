/**
 *  @file   larpandoracontent/LArCheating/CheatingCCLeptonClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CC_LEPTON_CLUSTER_CREATION_ALGORITHM_H
#define LAR_CHEATING_CC_LEPTON_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingCCLeptonClusterCreationAlgorithm class
 */
class CheatingCCLeptonClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingCCLeptonClusterCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCParticleToHitListMap;

    /**
     *  @brief  Create map between each (primary) MC particle and associated calo hits
     *
     *  @param  mcParticleToHitListMap to receive the mc particle to hit list map
     */
    void GetMCParticleToHitListMap(MCParticleToHitListMap &mcParticleToHitListMap) const;

    /**
     *  @brief  Simple mc particle collection, using main mc particle associated with each calo hit
     *
     *  @param  pCaloHit address of the calo hit
     *  @param  mcPrimaryMap the mapping between mc particles and their parents
     *  @param  mcParticleToHitListMap the mc particle to hit list map
     */
    void SimpleMCParticleCollection(const pandora::CaloHit *const pCaloHit, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap,
        MCParticleToHitListMap &mcParticleToHitListMap) const;

    /**
     *  @brief  Check whether mc particle is of a type specified for inclusion in cheated clustering
     *
     *  @param  pMCParticle the mc particle to hit list map
     *
     *  @return boolean
     */
    bool SelectMCParticlesForClustering(const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Create clusters based on information in the mc particle to hit list map
     *
     *  @param  mcParticleToHitListMap the mc particle to hit list map
     */
    void CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const;

    bool m_collapseToPrimaryMCParticles; ///< Whether to collapse mc particle hierarchies to primary particles
    std::string m_mcParticleListName;    ///< The mc particle list name, required if want to collapse mc particle hierarchy
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CC_LEPTON_CLUSTER_CREATION_ALGORITHM_H
