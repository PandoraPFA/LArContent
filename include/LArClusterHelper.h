/**
 *  @file   LArClusterHelper.h
 * 
 *  @brief  Header file for the cluster helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_HELPER_H
#define LAR_CLUSTER_HELPER_H 1

#include "Objects/Cluster.h"

namespace lar
{

/**
 *  @brief  ClusterQuality enum
 */
enum ClusterQuality
{
    METHOD_A = 0,  // NEED NAME
    METHOD_B,      // NEED NAME
    METHOD_C,      // NEED NAME
    METHOD_D       // NEED NAME
};

/**
 *  @brief  LArClusterHelper class
 */
class LArClusterHelper
{
public:

    /**
     *  @brief  Get length squared of cluster
     * 
     *  @param  pCluster address of the first cluster
     * 
     *  @return the length squared
     */
    static float GetLengthSquared( const pandora::Cluster* const pCluster );

    /**
     *  @brief  Get length of cluster
     * 
     *  @param  pCluster address of the first cluster
     * 
     *  @return the length
     */
    static float GetLength( const pandora::Cluster* const pCluster );

    /**
     *  @brief  Get number of layers spanned by cluster (1+Last-First)
     * 
     *  @param  pCluster address of the first cluster
     * 
     *  @return the layer span
     */

    static unsigned int GetLayerSpan( const pandora::Cluster* const pCluster );

    /**
     *  @brief  Fraction of occupied layers in cluster
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return float
     */
    static float GetLayerOccupancy(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in a pair of clusters
     * 
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     * 
     *  @return float
     */
    static float GetLayerOccupancy(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between the layer centroids of a pair of clusters
     * 
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     * 
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between a specified position vector and the layer centroids of a specified cluster
     * 
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     * 
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /** 
     * @brief Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    static void GetListOfCleanClusters(const ClusterQuality method, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector);

    /** 
     * @brief Populate cluster vector with clean clusters (method A)
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    static void GetListOfCleanClusters_MethodA(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector);

    /** 
     * @brief Populate cluster vector with clean clusters (method B)
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    static void GetListOfCleanClusters_MethodB(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector);

    /** 
     * @brief Populate cluster vector with clean clusters (method C)
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    static void GetListOfCleanClusters_MethodC(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector);

    /** 
     * @brief Populate cluster vector with clean clusters (method D)
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    static void GetListOfCleanClusters_MethodD(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector);

    /**
     *  @brief  Sort clusters by inner layer (then use SortByNOccupiedLayers method in event of a tie)
     * 
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByInnerLayer(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by number of occupied layers, and by inner layer, then energy in event of a tie
     * 
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByNOccupiedLayers(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by number of hits and by layer span, then energy in event of a tie
     * 
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByNHits(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Read the vertex helper settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_HELPER_H
