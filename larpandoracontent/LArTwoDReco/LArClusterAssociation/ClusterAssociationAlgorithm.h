/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h
 *
 *  @brief  Header file for the cluster association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_ASSOCIATION_ALGORITHM_H
#define LAR_CLUSTER_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ClusterAssociationAlgorithm class
 */
class ClusterAssociationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterAssociationAlgorithm();

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  ClusterAssociation class
     */
    class ClusterAssociation
    {
    public:
        pandora::ClusterSet m_forwardAssociations;  ///< The list of forward associations
        pandora::ClusterSet m_backwardAssociations; ///< The list of backward associations
    };

    typedef std::unordered_map<const pandora::Cluster *, ClusterAssociation> ClusterAssociationMap;

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
    virtual bool IsExtremalCluster(
        const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const = 0;

private:
    /**
     *  @brief  Unambiguous propagation
     *
     *  @param  pCluster address of the cluster to propagate
     *  @param  isForward whether propagation direction is forward
     *  @param  clusterAssociationMap the cluster association map
     */
    void UnambiguousPropagation(const pandora::Cluster *const pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Ambiguous propagation
     *
     *  @param  pCluster address of the cluster to propagate
     *  @param  isForward whether propagation direction is forward
     *  @param  clusterAssociationMap the cluster association map
     */
    void AmbiguousPropagation(const pandora::Cluster *const pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Update cluster association map to reflect an unambiguous cluster merge
     *
     *  @param  pClusterToEnlarge address of the cluster to be enlarged
     *  @param  pClusterToDelete address of the cluster to be deleted
     *  @param  isForwardMerge whether merge is forward (pClusterToEnlarge is forward-associated with pClusterToDelete)
     *  @param  clusterAssociationMap the cluster association map
     */
    void UpdateForUnambiguousMerge(const pandora::Cluster *const pClusterToEnlarge, const pandora::Cluster *const pClusterToDelete,
        const bool isForwardMerge, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Update cluster association map to reflect an ambiguous cluster merge
     *
     *  @param  pCluster address of the cluster to be cleared
     *  @param  clusterAssociationMap the cluster association map
     */
    void UpdateForAmbiguousMerge(const pandora::Cluster *const pCluster, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Navigate along cluster associations, from specified cluster, in specified direction
     *
     *  @param  clusterAssociationMap the cluster association map
     *  @param  pCluster address of cluster with which to begin search
     *  @param  isForward whether propagation direction is forward
     *  @param  pExtremalCluster to receive the extremal cluster
     *  @param  clusterSet to receive set of clusters traversed
     */
    void NavigateAlongAssociations(const ClusterAssociationMap &clusterAssociationMap, const pandora::Cluster *const pCluster,
        const bool isForward, const pandora::Cluster *&pExtremalCluster, pandora::ClusterSet &clusterSet) const;

    mutable bool m_mergeMade;

    bool m_resolveAmbiguousAssociations; ///< Whether to resolve ambiguous associations
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_ASSOCIATION_ALGORITHM_H
