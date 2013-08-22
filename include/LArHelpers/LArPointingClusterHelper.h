/**
 *  @file   LArContent/include/LArHelpers/LArClusterHelper.h
 * 
 *  @brief  Header file for the cluster helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_POINTING_CLUSTER_HELPER_H
#define LAR_POINTING_CLUSTER_HELPER_H 1

#include "Objects/Cluster.h"

#include "LArObjects/LArPointingCluster.h"

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
     *  @brief  Whether pointing vertex is adjacent to a given position
     * 
     *  @param  parentVertex the parent vertex position
     *  @param  daughtervertex the daughter pointing vertex
     * 
     *  @return boolean
     */
    static bool IsNode(const pandora::CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterVertex); 

    /**
     *  @brief  Whether pointing vertex is emitted from a given position
     * 
     *  @param  parentVertex the parent vertex position
     *  @param  daughtervertex the daughter pointing vertex
     * 
     *  @return boolean
     */
    static bool IsEmission(const pandora::CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterVertex);

    /**
     *  @brief  Get intersection of two vertices
     * 
     *  @param  firstVertex the first vertex
     *  @param  secondVertex the second vertex
     *  @param  intersectPosition the position of the intersection
     *  @param  firstDisplacement the displacement from first vertex to intersection
     *  @param  secondDisplacement the displacement from second vertex to intersection
     */
    static void GetIntersection(const LArPointingCluster::Vertex &firstVertex, const LArPointingCluster::Vertex &secondVertex,
        pandora::CartesianVector &intersectPosition, float &firstDisplacement, float &secondDisplacement);

    /**
     *  @brief  Get average direction of two vertices
     * 
     *  @param  firstVertex the first vertex
     *  @param  secondVertex the second vertex
     *  @param  averageDirection the average direction
     */
    static void GetAverageDirection(const LArPointingCluster::Vertex &firstVertex, const LArPointingCluster::Vertex &secondVertex,
        pandora::CartesianVector &averageDirection);

    /**
     *  @brief  Get intersection of two vertices
     * 
     *  @param  firstPosition the position of the first vertex
     *  @param  firstDirection the direction of the first vertex
     *  @param  secondPosition the position of the second vertex
     *  @param  secondDirection the direction of the second vertex
     *  @param  intersectPosition the position of the intersection
     *  @param  firstDisplacement the displacement from first vertex to intersection
     *  @param  secondDisplacement the displacement from second vertex to intersection
     */
    static void GetIntersection(const pandora::CartesianVector &firstPosition, const pandora::CartesianVector &firstDirection,
        const pandora::CartesianVector &secondPosition, const pandora::CartesianVector &secondDirection,
        pandora::CartesianVector &intersectPosition, float &firstDisplacement, float &secondDisplacement);

    /**
     *  @brief  Get intersection of vertex with target cluster
     * 
     *  @param  vertexCluster the vertex
     *  @param  pTargetCluster the target cluster
     *  @param  intersectPosition the position of the intersection
     *  @param  displacementL the longitudinal displacement 
     *  @param  displacementT the transverse displacement 
     */
    static void GetIntersection(const LArPointingCluster::Vertex &vertexCluster, const pandora::Cluster *const pTargetCluster,
        pandora::CartesianVector &intersectPosition, float &displacementL, float &displacementT);

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
