/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/IsolatedClusterMopUpAlgorithm.h
 * 
 *  @brief  Header file for the isolated cluster mop up algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_ISOLATED_CLUSTER_MOP_UP_ALGORITHM_H
#define LAR_ISOLATED_CLUSTER_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/ClusterMopUpBaseAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  IsolatedClusterMopUpAlgorithm class
 */
class IsolatedClusterMopUpAlgorithm : public ClusterMopUpBaseAlgorithm
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

    /**
     *  @brief  Default constructor
     */
    IsolatedClusterMopUpAlgorithm();

private:
    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters) const;

    /**
     *  @brief  Examine a list of clusters, identify and delete remnants; receive the list of newly available hits
     * 
     *  @param  clusterList the list of clusters to consider
     *  @param  caloHitList to receive the list of newly available hits
     */
    void DissolveClustersToHits(const pandora::ClusterList &clusterList, pandora::CaloHitList &caloHitList) const;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> CaloHitToClusterMap;

    /**
     *  @brief  Look for isolated hit additions, considering a list of candidate hits and a list of host clusters
     * 
     *  @param  caloHitList the list of hits to consider
     *  @param  clusterList the list of clusters to consider
     *  @param  caloHitToClusterMap to receive the calo hit to cluster map
     */
    void GetCaloHitToClusterMap(const pandora::CaloHitList &caloHitList, const pandora::ClusterList &clusterList, CaloHitToClusterMap &caloHitToClusterMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;

    unsigned int    m_maxCaloHitsInCluster;     ///< The maximum number of hits in a cluster to be dissolved
    float           m_maxHitClusterDistance;    ///< The maximum hit to cluster distance for isolated hit merging
    bool            m_addHitsAsIsolated;        ///< Whether to add hits to clusters as "isolated" (don't contribute to spatial properties)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *IsolatedClusterMopUpAlgorithm::Factory::CreateAlgorithm() const
{
    return new IsolatedClusterMopUpAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_ISOLATED_CLUSTER_MOP_UP_ALGORITHM_H
