/**
 *  @file   LArContent/include/LArClusterAssociation/ClusterMergingAlgorithm.h
 * 
 *  @brief  Header file for the cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
#define LAR_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArOverlapTensor.h"

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

    typedef std::map<pandora::Cluster*, bool> ClusterVetoMap;
    typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

    typedef OverlapTensor<bool> ClusterAssociationMatrix;

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
     *  @param  clusterMergingMatrix the matrix of cluster associations
     */
    virtual void FillAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const = 0;

  
private:

    /**
     *  @brief  Determine whether two clusters are associated and should be merged
     * 
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     *
     *  @return boolean
     */

    bool AreClustersAssociated(pandora::Cluster *pCluster1, pandora::Cluster *pCluster2, const ClusterAssociationMatrix &clusterAssociationMatrix);

   /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     * 
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  pCurrentCluster pointer to the current cluster 
     *  @param  clusterMergeMap the map of clusters to be merged
     *  @param  clusterVetoMap the map of clusters that have already been merged
     *  @param  associatedClusterList the list of associated clusters
     */
    void CollectAssociatedClusters(pandora::Cluster *pSeedCluster, pandora::Cluster *pCurrentCluster, const ClusterMergeMap &clusterMergeMap,
        const ClusterVetoMap &clusterVetoMap, pandora::ClusterList& associatedClusterList) const;
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
