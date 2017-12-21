/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayShowerMatchingAlg.h
 * 
 *  @brief  Header file for the cosmic ray shower matching cheater class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_COSMIC_RAY_SHOWER_MATCHING_ALG_H
#define LAR_CHEATING_COSMIC_RAY_SHOWER_MATCHING_ALG_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingCosmicRayShowerMatchingAlg class
 */
class CheatingCosmicRayShowerMatchingAlg : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the list of candidate clusters for matching with existing pfos
     *
     *  @param  candidateClusterList to receive the list of candidate clusters
     */
    void GetCandidateClusters(pandora::ClusterList &candidateClusterList) const;

    /**
     *  @brief  Perform cosmic ray shower matching for a specific cluster in a pfo
     *
     *  @param  pPfo the pfo of interest
     *  @param  pPfoCluster the pfo cluster of interest
     *  @param  candidateClusterList the list of candidate clusters
     */
    void CosmicRayShowerMatching(const pandora::ParticleFlowObject *const pPfo, const pandora::Cluster *const pPfoCluster,
        const pandora::ClusterList &candidateClusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_inputPfoListName;           ///< The input pfo list name
    pandora::StringVector   m_inputClusterListNames;      ///< The input cluster list names
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_COSMIC_RAY_SHOWER_MATCHING_ALG_H
