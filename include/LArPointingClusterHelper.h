/**
 *  @file   LArClusterHelper.h
 * 
 *  @brief  Header file for the cluster helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_POINTING_CLUSTER_HELPER_H
#define LAR_POINTING_CLUSTER_HELPER_H 1

#include "Objects/Cluster.h"

#include "LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  LArPointingClusterHelper class
 */
class LArPointingClusterHelper
{
public:

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
     *  @brief  Get intersection of two vertices
     * 
     *  @param  firstVertex the first vertex
     *  @param  secondVertex the second vertex
     *  @param  intersectPosition the position of the intersection
     *  @param  isPhysical do they obey causality
     */
    static void GetIntersection( const LArPointingCluster::Vertex& firstVertex, const LArPointingCluster::Vertex& secondVertex, pandora::CartesianVector& intersectPosition, bool& isPhysical );

    /**
     *  @brief  Get average direction of two vertices
     * 
     *  @param  firstVertex the first vertex
     *  @param  secondVertex the second vertex
     *  @param  averageDirection the average direction
     */
    static void GetAverageDirection( const LArPointingCluster::Vertex& firstVertex, const LArPointingCluster::Vertex& secondVertex, pandora::CartesianVector& averageDirection );

    /**
     *  @brief  Get intersection of two vertices
     * 
     *  @param  firstPosition the position of the first vertex
     *  @param  firstDirection the direction of the first vertex
     *  @param  secondPosition the position of the second vertex
     *  @param  secondDirection the direction of the second vertex
     *  @param  intersectPosition the position of the intersection
     *  @param  isPhysical do they obey causality
     */
    static void GetIntersection( const pandora::CartesianVector& firstPosition, const pandora::CartesianVector& firstDirection,
        const pandora::CartesianVector& secondPosition, const pandora::CartesianVector& secondDirection,
        pandora::CartesianVector& intersectPosition, bool& isPhysical );

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

#endif // #ifndef LAR_POINTING_CLUSTER_HELPER_H
