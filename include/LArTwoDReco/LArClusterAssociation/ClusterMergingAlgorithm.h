/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the cluster merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
#define LAR_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ClusterMergingAlgorithm class
 */
class ClusterMergingAlgorithm : public pandora::Algorithm
{
protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const pandora::Cluster*, bool> ClusterVetoMap;
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    virtual void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const = 0;

    /**
     *  @brief  Form associations between pointing clusters
     *
     *  @param  clusterVector the vector of clean clusters
     *  @param  clusterMergeMap the matrix of cluster associations
     */
    virtual void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const = 0;

    /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     *
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  clusterMergeMap the map of cluster associations
     *  @param  associatedClusterList the output list of associated clusters
     */
    void CollectAssociatedClusters(pandora::Cluster *pSeedCluster, const ClusterMergeMap &clusterMergeMap, pandora::ClusterList& associatedClusterList) const;

   /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     *
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  pCurrentCluster pointer to the current cluster
     *  @param  clusterMergeMap the map of cluster associations
     *  @param  clusterVetoMap the map of clusters that have already been merged
     *  @param  associatedClusterList the output list of associated clusters
     */
    void CollectAssociatedClusters(pandora::Cluster *pSeedCluster, pandora::Cluster *pCurrentCluster, const ClusterMergeMap &clusterMergeMap,
        const ClusterVetoMap &clusterVetoMap, pandora::ClusterList& associatedClusterList) const;

    std::string     m_inputClusterListName;     ///< The name of the input cluster list. If not specified, will access current list.
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
