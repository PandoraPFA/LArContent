/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster merging algorithm.
 *
 *  $Log: $
 */

#ifndef LAR_CHEATING_CLUSTER_MERGING_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_MERGING__ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingClusterMergingAlgorithm class
 */
class CheatingClusterMergingAlgorithm : public pandora::Algorithm
{
public:
/**
 *  @brief  Default constructor
 */
    CheatingClusterMergingAlgorithm();
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Cheated Cluster Merging. Use MC to match clusters based on the main MC particle.
     *
     *  @param  pClusterList the list of clusters
     *  @param  listName the name of the current cluster list (ClustersU, ClustersV, ClustersW normally)
     */
    void CheatedClusterMerging(const pandora::ClusterList *const pClusterList, const std::string &listName) const;

     /**
      *  @brief  Get the MC particle for a given cluster, caching to a map.
      *
      *  @param  cluster         The current cluster to lookup
      *  @param  clusterToMCMap  Map from Cluster to MC to cache results.
      */
    const pandora::MCParticle* GetMCForCluster(const pandora::Cluster *const cluster, std::map<const pandora::Cluster*,
        const pandora::MCParticle*> &clusterToMCMap) const;

    /**
     *  @brief  If a cluster is valid to use: Is a shower tagged cluster, and not been used yet.
     *
     *  @param  cluster         The current cluster to lookup
     *  @param  clusterIsUsed  Map from Cluster to bool to check if a cluster has been used yet.
     */
    bool IsValidToUse(const pandora::Cluster *const cluster, std::map<const pandora::Cluster*, bool> &clusterIsUsed) const;

    pandora::StringVector  m_inputClusterListNames; ///< The names of the input cluster lists.
    std::string            m_mcParticleListName;    ///< Input MC particle list name.
    float                  m_minNCaloHits;          ///< The minimum number of hits for a cluster to be deemed true for IsAvailableToUse.
    

};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CLUSTER_MERGING_ALGORITHM_H
