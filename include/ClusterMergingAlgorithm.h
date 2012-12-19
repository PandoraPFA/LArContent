/**
 *  @file   ClusterMergingAlgorithm.h
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     * 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters( const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector );
    
    /**
     *  @brief  Determine whether two clusters are associated and should be merged
     * 
     *  @param  pClusterI the address of the first cluster
     *  @param  pClusterJ the address of the second cluster
     *
     *  @return boolean
     */
    bool AreClustersAssociated( const pandora::Cluster* const pClusterI, const pandora::Cluster* const pClusterJ );

    /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     * 
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  pCurrentCluster pointer to the current cluster 
     *  @param  associatedClusterList the list of associated clusters
     */
    void CollectAssociatedClusters( pandora::Cluster* pSeedCluster, pandora::Cluster* pCurrentCluster, pandora::ClusterList& associatedClusterList );


    typedef std::map<pandora::Cluster*, bool> ClusterVetoMap;
    ClusterVetoMap m_clusterVetoMap;

    typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;
    ClusterMergeMap m_clusterMergeMap;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
