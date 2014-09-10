/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h
 *
 *  @brief  Header file for the cluster growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_GROWING_ALGORITHM_H
#define LAR_CLUSTER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterGrowingAlgorithm class
 */
class ClusterGrowingAlgorithm : public pandora::Algorithm
{
protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

    /**
     *  @brief  Populate cluster vector with the subset of clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  cleanClusters the output vector of clean clusters
     */
    virtual void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &cleanClusters) const = 0;

    /**
     *  @brief Select seed clusters for growing
     *
     *  @param cleanClusters the input vector of clean clusters
     *  @param seedClusters the output vector of seed clusters
     */
    virtual void GetListOfSeedClusters(const pandora::ClusterVector &cleanClusters, pandora::ClusterVector &seedClusters) const = 0;

private:

    /**
     *  @brief Get List of non-seed clusters
     *
     *  @param inputClusters the input vector of clean clusters
     *  @param seedClusters the input vector of seed clusters
     *  @param nonSeedClusters the output vector of non-seed clusters
     */
    void GetListOfNonSeedClusters(const pandora::ClusterVector &inputClusters, const pandora::ClusterVector &seedClusters,
        pandora::ClusterVector &nonSeedClusters) const;

    /**
     *  @brief Identify a set of cluster merges
     *
     *  @param seedClusters the input vector of seed clusters
     *  @param nonSeedClusters the input vector of non-seed clusters
     *  @param clusterMergeMap the output map of cluster merges
     */
    void PopulateClusterMergeMap(const pandora::ClusterVector &seedClusters, const pandora::ClusterVector &nonSeedClusters,
        ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief Merge clusters
     *
     *  @param clusterMergeMap the map of clusters to be merged
     */
    void MergeClusters(const ClusterMergeMap &clusterMergeMap) const;

    std::string     m_inputClusterListName;  ///< The name of the input cluster list. If not specified, will access current list.

    float           m_maxClusterSeparation;  ///< Maximum distance at which clusters can be joined
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_GROWING_ALGORITHM_H
