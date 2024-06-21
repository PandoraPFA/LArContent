/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/BoundedRegionClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the bounded region cluster merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_BOUNDED_REGION_CLUSTER_MERGING_ALGORITHM_H
#define LAR_BOUNDED_REGION_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  BoundedRegionClusterMergingAlgorithm class
 */
class BoundedRegionClusterMergingAlgorithm : public ClusterMergingAlgorithm
{
public:
    BoundedRegionClusterMergingAlgorithm();

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::Cluster *, pandora::ClusterList> ClusterMergeMap;

    /**
     *  @brief Populate cluster vector with all clusters starting in the defined region
     *
     *  @param  pClusterList pointer to the list of all 2D clusters
     *  @param  clusterVector to receive the clusters within the bounded region
     */
    void GetListOfBoundedRegionClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief Compare two clusters and decide if they should be associated with each other
     *
     *  @param pCluster1 pointer to the first cluster
     *  @param pCluster2 pointer to the second cluster
     *
     *  @return whether the clusters should be associated or not
     */
    bool AreClustersAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Form associations between pointing clusters
     *
     *  @param  clusterVector the vector of clean clusters
     *  @param  clusterMergeMap the matrix of cluster associations
     */
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

private:
    float m_xRegionMin;            ///< Minimum x value of the bounded region box
    float m_xRegionMax;            ///< Maximum x value of the bounded region box
    float m_zRegionMin;            ///< Minimum z value (wire dimension) of the bounded region box
    float m_zRegionMax;            ///< Maximum z value (wire dimension) of the bounded region box
    float m_maxDistance;           ///< Maximum distance below which clusters can be associated
    unsigned int m_minClusterHits; ///< Threshold on the size of clusters to be considered for merging
};

} // namespace lar_content

#endif // #ifndef LAR_BOUNDED_REGION_CLUSTER_MERGING_ALGORITHM_H
