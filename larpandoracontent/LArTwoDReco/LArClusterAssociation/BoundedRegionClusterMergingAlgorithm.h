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

    void GetListOfBoundedRegionClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    bool AreClustersAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const {};

    /**
     *  @brief  Form associations between pointing clusters
     *
     *  @param  clusterVector the vector of clean clusters
     *  @param  clusterMergeMap the matrix of cluster associations
     */
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Merge associated clusters
     *
     *  @param  clusterVector the vector of clean clusters
     *  @param  clusterMergeMap the matrix of cluster associations
     */
//    void MergeClusters(pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     *
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  clusterMergeMap the map of cluster associations
     *  @param  associatedClusterList the output list of associated clusters
     */
//    void CollectAssociatedClusters(const pandora::Cluster *const pSeedCluster, const ClusterMergeMap &clusterMergeMap,
//        pandora::ClusterList &associatedClusterList) const;

    /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     *
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  pCurrentCluster pointer to the current cluster
     *  @param  clusterMergeMap the map of cluster associations
     *  @param  clusterVetoList the list of clusters that have already been merged
     *  @param  associatedClusterList the output list of associated clusters
     */
//    void CollectAssociatedClusters(const pandora::Cluster *const pSeedCluster, const pandora::Cluster *const pCurrentCluster,
//        const ClusterMergeMap &clusterMergeMap, const pandora::ClusterSet &clusterVetoList, pandora::ClusterList &associatedClusterList) const;

    /**
     *  @brief  Sort the selected clusters, so that they have a well-defined ordering
     *
     *  @param  inputClusters the input vector of clusters
     *  @param  outputClusters the output vector of clusters
     */
//    void GetSortedListOfCleanClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &outputClusters) const;

private:

    float m_xRegionMin;
    float m_xRegionMax;
    float m_zRegionMin;
    float m_zRegionMax;
    float m_maxDistance;
    unsigned int m_minClusterHits;
};

} // namespace lar_content

#endif // #ifndef LAR_BOUNDED_REGION_CLUSTER_MERGING_ALGORITHM_H
