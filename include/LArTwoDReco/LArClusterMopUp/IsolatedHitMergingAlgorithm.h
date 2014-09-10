/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h
 * 
 *  @brief  Header file for the isolated hit merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_ISOLATED_HIT_MERGING_ALGORITHM_H
#define LAR_ISOLATED_HIT_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  IsolatedHitMergingAlgorithm class
 */
class IsolatedHitMergingAlgorithm : public ClusterMopUpAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters, const ClusterToListNameMap &clusterToListNameMap) const;

    /**
     *  @brief  Examine a list of clusters, identify and delete remnants; receive the list of newly available hits
     * 
     *  @param  clusterList the list of clusters to consider
     *  @param  clusterToListNameMap the cluster to list name map
     *  @param  caloHitList to receive the list of newly available hits
     */
    void DissolveClustersToHits(const pandora::ClusterList &clusterList, const ClusterToListNameMap &clusterToListNameMap, pandora::CaloHitList &caloHitList) const;

    typedef std::map<pandora::CaloHit*, pandora::Cluster*> CaloHitToClusterMap;

    /**
     *  @brief  Look for isolated hit additions, considering a list of candidate hits and a list of host clusters
     * 
     *  @param  caloHitList the list of hits to consider
     *  @param  clusterList the list of clusters to consider
     *  @param  caloHitToClusterMap to receive the calo hit to cluster map
     */
    void GetCaloHitToClusterMap(const pandora::CaloHitList &caloHitList, const pandora::ClusterList &clusterList, CaloHitToClusterMap &caloHitToClusterMap) const;

    /**
     *  @brief  Get closest distance between a specified calo hit and a non-isolated hit in a specified cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  pCaloHit address of the calo hit
     * 
     *  @return The closest distance between the calo hit and a non-isolated hit in the cluster
     */
    float GetDistanceToHit(const pandora::Cluster *const pCluster, const pandora::CaloHit *const pCaloHit) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_maxCaloHitsInCluster;     ///< The maximum number of hits in a cluster to be dissolved
    int             m_hitLayerSearchWindow;     ///< The layer search window used (for speed) when calculating hit to cluster distances
    float           m_maxHitClusterDistance;    ///< The maximum hit to cluster distance for isolated hit merging
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *IsolatedHitMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new IsolatedHitMergingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_ISOLATED_HIT_MERGING_ALGORITHM_H
