/**
 *  @file   LArContent/include/LArHelpers/LArClusterHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_HELPER_H
#define LAR_CLUSTER_HELPER_H 1

#include "Objects/Cluster.h"

namespace lar_content
{

/**
 *  @brief  LArClusterHelper class
 */
class LArClusterHelper
{
public:
    /**
     *  @brief  Get the hit type associated with a two dimensional cluster
     * 
     *  @param  pCluster the address of the cluster
     * 
     *  @return the cluster hit type
     */
    static pandora::HitType GetClusterHitType(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get length squared of cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the length squared
     */
    static float GetLengthSquared(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get length of cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the length
     */
    static float GetLength(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get energy of cluster, based on length
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the energy
     */
    static float GetEnergyFromLength(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get number of layers spanned by cluster (1+Last-First)
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the layer span
     */
    static unsigned int GetLayerSpan(const pandora::Cluster *const pCluster);

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
     *  @brief  Get closest distance between clusters in a pair of cluster lists
     *
     *  @param  clusterList1 the first cluster list
     *  @param  clusterList2 the second cluster list
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2);

    /**
     *  @brief  Get closest distance between a specified cluster and list of clusters
     *
     *  @param  pCluster address of the input cluster
     *  @param  clusterList list of input clusters
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster, const pandora::ClusterList &clusterList);

    /**
     *  @brief  Get closest distance between a pair of clusters
     *
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between a specified position vector and the hits in a specified cluster
     *
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get closest position on a cluster to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get positions of the two most distant calo hits in a 2D cluster (ordered by Z)
     *
     *  @param  pCluster address of the cluster
     *  @param  the inner extremal position
     *  @param  the outer extremal position
     */
    static void GetExtremalCoordinatesXZ(const pandora::Cluster *const pCluster, pandora::CartesianVector &innerCoordinate, 
        pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief  Get minimum and maximum X, Y and Z positions of the calo hits in a 2D cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  the minimum positions (x,y,z)
     *  @param  the maximum positions (x,y,z)
     */
    static void GetClusterSpanXZ(const pandora::Cluster *const pCluster, pandora::CartesianVector &minimumCoordinate, 
        pandora::CartesianVector &maximumCoordinate);

    /**
     *  @brief  Get minimum and maximum x of the calo hits in a 2D cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  the minimum position of x
     *  @param  the maximum position of x
     */
    static void GetClusterSpanX(const pandora::Cluster *const pCluster, float &xmin, float &xmax);

    /**
     *  @brief  Get upper and lower Z positions of the calo hits in a 2D cluster in range xmin to xmax
     *
     *  @param  pCluster address of the cluster
     *  @param  xmin for range in x
     *  @param  xmax for range in x
     *  @param  zmin the lower z for this range of x
     *  @param  zmax the upper z for this range in x
     */
    static void GetClusterSpanZ(const pandora::Cluster *const pCluster, const float xmin, const float xmax, float &zmin, float &zmax);

    /**
     *  @brief  Get average Z positions of the calo hits in a 2D cluster in range xmin to xmax
     *
     *  @param  pCluster address of the cluster
     *  @param  xmin for range in x
     *  @param  xmax for range in x
     */
    static float GetAverageZ(const pandora::Cluster *const pCluster, const float xmin, const float xmax);

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
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_HELPER_H
