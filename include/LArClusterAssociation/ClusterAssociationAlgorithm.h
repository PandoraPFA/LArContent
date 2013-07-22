/**
 *  @file   LArContent/include/LArClusterAssociation/ClusterAssociationAlgorithm.h
 * 
 *  @brief  Header file for the cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_ASSOCIATION_ALGORITHM_H
#define LAR_CLUSTER_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Helpers/ClusterHelper.h"

namespace lar
{

/**
 *  @brief  ClusterAssociationAlgorithm class
 */
class ClusterAssociationAlgorithm : public pandora::Algorithm
{
protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  ClusterAssociation class
     */
    class ClusterAssociation
    {
    public:
        pandora::ClusterList    m_forwardAssociations;      ///< The list of forward associations
        pandora::ClusterList    m_backwardAssociations;     ///< The list of backward associations
    };

    typedef std::map<pandora::Cluster*, ClusterAssociation> ClusterAssociationMap;

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     * 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    virtual void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const = 0;

    /**
     *  @brief  Populate the cluster association map
     * 
     *  @param  clusterVector the cluster vector
     *  @param  clusterAssociationMap to receive the populated cluster association map
     */
    virtual void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const = 0;

    /**
     *  @brief  Determine which of two clusters is extremal
     * 
     *  @param  isForward whether propagation direction is forward
     *  @param  pCurrentCluster current extremal cluster
     *  @param  pTestCluster potential extremal cluster
     * 
     *  @return boolean
     */
    virtual bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const = 0;

private:
    /**
     *  @brief  Unambiguous propagation
     * 
     *  @param  pCluster address of the cluster to propagate
     *  @param  isForward whether propagation direction is forward
     *  @param  clusterAssociationMap the cluster association map
     */
    void UnambiguousPropagation(pandora::Cluster *pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Ambiguous propagation
     * 
     *  @param  pCluster address of the cluster to propagate
     *  @param  isForward whether propagation direction is forward
     *  @param  clusterAssociationMap the cluster association map
     */
    void AmbiguousPropagation(pandora::Cluster *pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Update cluster association map to reflect an unambiguous cluster merge
     * 
     *  @param  pClusterToEnlarge address of the cluster to be enlarged
     *  @param  pClusterToDelete address of the cluster to be deleted
     *  @param  isForwardMerge whether merge is forward (pClusterToEnlarge is forward-associated with pClusterToDelete)
     *  @param  clusterAssociationMap the cluster association map
     */
    void UpdateForUnambiguousMerge(pandora::Cluster *pClusterToEnlarge, pandora::Cluster *pClusterToDelete, const bool isForwardMerge,
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Update cluster association map to reflect an ambiguous cluster merge
     * 
     *  @param  pClusterToEnlarge address of the cluster to be enlarged
     *  @param  pClusterToDelete address of the cluster to be deleted
     *  @param  isForwardMerge whether merge is forward (pClusterToEnlarge is forward-associated with pClusterToDelete)
     *  @param  clusterAssociationMap the cluster association map
     */
    void UpdateForAmbiguousMerge(pandora::Cluster *pClusterToEnlarge, pandora::Cluster *pClusterToDelete, const bool isForwardMerge,
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Navigate along cluster associations, from specified cluster, in specified direction
     * 
     *  @param  clusterAssociationMap the cluster association map
     *  @param  pCluster address of cluster with which to begin search
     *  @param  isForward whether propagation direction is forward
     *  @param  pExtremalCluster to receive the extremal cluster
     *  @param  clusterList to receive list of clusters traversed
     */
    void NavigateAlongAssociations(const ClusterAssociationMap &clusterAssociationMap, pandora::Cluster *pCluster, const bool isForward,
        pandora::Cluster *&pExtremalCluster, pandora::ClusterList &clusterList) const;

    mutable bool m_mergeMade;

    bool         m_resolveAmbiguousAssociations;        ///< Whether to resolve ambiguous associations
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_ASSOCIATION_ALGORITHM_H
