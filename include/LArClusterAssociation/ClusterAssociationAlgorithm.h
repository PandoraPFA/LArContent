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
    virtual void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Populate the cluster association map
     * 
     *  @param  clusterVector the cluster vector
     *  @param  clusterAssociationMap to receive the populated cluster association map
     */
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Determine whether two clusters are associated
     * 
     *  @param  pInnerCluster address of the inner cluster
     *  @param  pOuterCluster address of the outer cluster
     * 
     *  @return whether the clusters are associated
     */
    virtual bool AreClustersAssociated(const pandora::Cluster *const pInnerCluster, const pandora::Cluster *const pOuterCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     * 
     *  @param  innerClusterEnd inner cluster end position
     *  @param  outerClusterStart outer cluster start position
     *  @param  hitSizeX hit size x
     *  @param  hitSizeZ hit size z
     *  @param  innerFit inner cluster fit result
     *  @param  outerFit outer cluster fit result
     * 
     *  @return whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::CartesianVector &innerClusterEnd, const pandora::CartesianVector &outerClusterStart, const float hitSizeX,
        const float hitSizeZ, const pandora::ClusterHelper::ClusterFitResult &innerFit, const pandora::ClusterHelper::ClusterFitResult &outerFit) const;

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

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterAssociationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterAssociationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_ASSOCIATION_ALGORITHM_H
