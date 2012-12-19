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

#include "LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  LArClusterHelper class
 */
class LArClusterHelper
{
public:
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
     *  @brief  Whether specified parent and daughter vertices form a node
     * 
     *  @param  parentVertex the position of the parent cluster vertex
     *  @param  daughterVertex the position of the daughter cluster vertex
     * 
     *  @return boolean
     */
    static bool IsNode(const pandora::CartesianVector &parentVertex, const pandora::CartesianVector &daughterVertex);

    /**
     *  @brief  Whether parent and daughter vertices represent an emission from the parent
     * 
     *  @param  parentPointingVertex the parent pointing vertex
     *  @param  daughterVertex the position of the daughter vertex
     * 
     *  @return boolean
     */
    static bool IsEmission(const LArPointingCluster::Vertex &parentPointingVertex, const pandora::CartesianVector &daughterVertex);

    /**
     *  @brief  Whether parent and daughter vertices represent an emission from the parent
     * 
     *  @param  parentVertex the position of the parent vertex
     *  @param  daughterPointingVertex the daughter pointing vertex
     * 
     *  @return boolean
     */
    static bool IsEmitted(const pandora::CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterPointingVertex);

    /**
     *  @brief  Whether vertex and pointing vertex are compatible with an emission
     * 
     *  @param  vertex the position of the vertex
     *  @param  pointingVertex the pointing vertex
     * 
     *  @return boolean
     */
    static bool IsPointing(const pandora::CartesianVector &vertex, const LArPointingCluster::Vertex &pointingVertex);

    /**
     *  @brief  Read the vertex helper settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    static float    m_maxNodeRadiusSquared;                 ///< 
    static float    m_maxPointingLongitudinalDistance;      ///< 
    static float    m_minPointingLongitudinalDistance;      ///<
    static float    m_maxPointingTransverseDistance;        ///< 
    static float    m_pointingAngularAllowance;             ///< 
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_HELPER_H
