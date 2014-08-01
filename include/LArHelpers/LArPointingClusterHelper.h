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
     *  @brief  Calculate distance squared between inner and outer vertices of pointing cluster
     * 
     *  @param  pointingCluster  the input pointing cluster
     * 
     *  @return float  the distance squared
     */
    static float GetLengthSquared(const LArPointingCluster &pointingCluster);

    /**
     *  @brief  Calculate distance squared between inner and outer vertices of pointing cluster
     * 
     *  @param  pointingCluster  the input pointing cluster
     * 
     *  @return float  the distance
     */
    static float GetLength(const LArPointingCluster &pointingCluster);

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
     *  @brief  Get projected distance between a cluster and a pointing cluster vertex
     * 
     *  @param  pointingVertex the pointing vertex
     *  @param  pCluster address of the cluster
     * 
     *  @return the projected distance
     */
    static float GetProjectedDistance(const LArPointingCluster::Vertex &pointingVertex, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get projected position on a cluster from a pointing cluster vertex
     * 
     *  @param  pointingVertex the pointing vertex
     *  @param  pCluster address of the cluster
     * 
     *  @return the projected position
     */
    static pandora::CartesianVector GetProjectedPosition(const LArPointingCluster::Vertex &pointingVertex, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get projected position on a cluster from a specified position and direction
     * 
     *  @param  initialPosition the initial position of the cluster
     *  @param  initialDirection the initial direction of the cluster
     *  @param  pCluster address of the cluster
     * 
     *  @return the projected position
     */
    static pandora::CartesianVector GetProjectedPosition(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Given a pair of pointing clusters, receive the closest pair of vertices (one vertex coming from each cluster)
     * 
     *  @param  pointingClusterI the first pointing cluster
     *  @param  pointingClusterJ the second pointing cluster
     *  @param  closestVertexI to receive the relevant vertex from the first pointing cluster
     *  @param  closestVertexJ to receive the relevant vertex from the second pointing cluster
     */
    static void GetClosestVertices(const LArPointingCluster &pointingClusterI, const LArPointingCluster &pointingClusterJ,
        LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ);

    /**
     *  @brief  Given a pair of pointing clusters, receive the pair of vertices with smallest x-separation (one vertex coming from each cluster)
     * 
     *  @param  pointingClusterI the first pointing cluster
     *  @param  pointingClusterJ the second pointing cluster
     *  @param  closestVertexI to receive the relevant vertex from the first pointing cluster
     *  @param  closestVertexJ to receive the relevant vertex from the second pointing cluster
     */
    static void GetClosestVerticesInX(const LArPointingCluster &pointingClusterI, const LArPointingCluster &pointingClusterJ,
        LArPointingCluster::Vertex &closestVertexI, LArPointingCluster::Vertex &closestVertexJ);

    /**
     *  @brief  Calculate impact parameters between a pointing cluster vertex and a target position
     * 
     *  @param  pointingVertex the pointing vertex
     *  @param  targetPosition the target position
     *  @param  longitidinal the longitudinal displacement
     *  @param  transverse the transverse displacement
     */
    static void GetImpactParameters(const LArPointingCluster::Vertex &pointingVertex, const LArPointingCluster::Vertex &targetVertex,
        float &longitudinal, float &transverse);

    /**
     *  @brief  Calculate impact parameters between a pointing cluster vertex and a target position
     * 
     *  @param  pointingVertex the pointing vertex
     *  @param  targetPosition the target position
     *  @param  longitidinal the longitudinal displacement
     *  @param  transverse the transverse displacement
     */
    static void GetImpactParameters(const LArPointingCluster::Vertex &pointingVertex, const pandora::CartesianVector &targetPosition, 
        float &longitudinal, float &transverse);

    /**
     *  @brief  Calculate impact parameters of a specified position and direction to a target position
     * 
     *  @param  initialPosition the initial position of the cluster
     *  @param  initialDirection the initial direction of the cluster
     *  @param  targetPosition the target position
     *  @param  longitidinal the longitudinal displacement
     *  @param  transverse the transverse displacement
     */
    static void GetImpactParameters(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
        const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse);

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

    /**
     *  @brief  Simple and fast vertex selection, choosing best vertex from a specified list to represent a set of pointing clusters
     * 
     *  @param  vertexList the candidate vertex list
     *  @param  pointingClusterList the pointing cluster list
     * 
     *  @return the best vertex estimate
     */
    static LArPointingCluster::Vertex GetBestVertexEstimate(const LArPointingClusterVertexList &vertexList, const LArPointingClusterList &pointingClusterList);

private:
    /**
     *  @brief  Collect cluster vertices, from a provided input list, associated with a specified vertex
     * 
     *  @param  vertex the vertex
     *  @param  inputList the input list of pointing clusters
     *  @param  outputList to receive the output list of cluster vertices associated with the specified vertex
     */
    static void CollectAssociatedClusters(const LArPointingCluster::Vertex &vertex, const LArPointingClusterList &inputList, LArPointingClusterVertexList &outputList);

    /**
     *  @brief  Get an estimate of the energy associated with a specified vertex
     * 
     *  @param  vertex the vertex
     *  @param  clusterVertices the list of cluster vertices associated with the specified vertex
     * 
     *  @return the energy associated with a specified vertex
     */
    static float GetAssociatedEnergy(const LArPointingCluster::Vertex &vertex, const LArPointingClusterVertexList &clusterVertices);

    static float    m_maxNodeRadiusSquared;                 ///< 
    static float    m_maxPointingLongitudinalDistance;      ///< 
    static float    m_minPointingLongitudinalDistance;      ///<
    static float    m_maxPointingTransverseDistance;        ///< 
    static float    m_pointingAngularAllowance;             ///< 
    static float    m_projectionAngularAllowance;           ///< 
};

} // namespace lar

#endif // #ifndef LAR_POINTING_CLUSTER_HELPER_H
