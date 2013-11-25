/**
 *  @file   LArContent/include/LArHelpers/LArVertexHelper.h
 * 
 *  @brief  Header file for the vertex helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_HELPER_H
#define LAR_VERTEX_HELPER_H 1

#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"

#include <map>
#include <string>

namespace lar
{

/**
 *  @brief  ClusterDirection enum
 */
enum ClusterDirection
{
    FORWARD = 0,
    BACKWARD,
    DIRECTION_AMBIGUOUS,
    DIRECTION_UNKNOWN
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArVertexHelper class
 */
class LArVertexHelper
{
public:
    /**
     *  @brief  Add a named vertex to the name to vertex map
     * 
     *  @param  vertexName the name/label for the vertex
     *  @param  vertex the vertex
     */
    static void AddVertex(const std::string &vertexName, const pandora::CartesianVector &vertex);

    /**
     *  @brief  Remove a named vertex from the name to vertex map
     * 
     *  @param  vertexName the name/label for the vertex
     */
    static void RemoveVertex(const std::string &vertexName);

    /**
     *  @brief  Update a named vertex in the name to vertex map
     * 
     *  @param  vertexName the name/label for the vertex
     *  @param  vertex the updated vertex
     */
    static void UpdateVertex(const std::string &vertexName, const pandora::CartesianVector &vertex);

    /**
     *  @brief  Set a named vertex to be the current vertex
     * 
     *  @param  vertexName the name/label for the current vertex
     */
    static void SetCurrentVertex(const std::string &vertexName);

    /**
     *  @brief  Get a named vertex from the name to vertex map
     * 
     *  @param  vertexName the name/label for the vertex
     * 
     *  @return the named vertex
     */
    static const pandora::CartesianVector &GetVertex(const std::string &vertexName);

    /**
     *  @brief  Get the current vertex, if set
     * 
     *  @return the current vertex
     */
    static const pandora::CartesianVector &GetCurrentVertex(); 

    /**
     *  @brief  Get the current vertex name, if set
     * 
     *  @return the current vertex name
     */
    static const std::string GetCurrentVertexName();

    /**
     *  @brief  Ask if the current vertex exists
     *
     *  @return boolean
     */
    static bool DoesCurrentVertexExist();

    /**
     *  @brief  Ask if a particular vertex exists
     *
     *  @return boolean
     */
    static bool DoesVertexExist(const std::string &vertexName);

    /**
     *  @brief  Ask for cluster propagation direction
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return the cluster direction
     */
    static ClusterDirection GetDirectionInZ(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Ask if cluster is propagation forwards in Z
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static bool IsForwardInZ(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Ask if cluster is propagation forwards in Z, wrt a specified vertex
     * 
     *  @param  pCluster address of the cluster
     *  @param  specifiedVertex the specified vertex
     *
     *  @return boolean
     */
    static bool IsForwardInZ(const pandora::Cluster *const pCluster, const pandora::CartesianVector &specifiedVertex);

    /**
     *  @brief  Ask if cluster is propagation backward in Z
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static bool IsBackwardInZ(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Ask if cluster is propagation backward in Z, wrt a specified vertex
     * 
     *  @param  pCluster address of the cluster
     *  @param  specifiedVertex the specified vertex
     *
     *  @return boolean
     */
    static bool IsBackwardInZ(const pandora::Cluster *const pCluster, const pandora::CartesianVector &specifiedVertex);

    /**
     *  @brief  Ask if cluster is propagation forwards in Z
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static bool IsForwardInZ3D(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Ask if cluster is propagation backward in Z
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static bool IsBackwardInZ3D(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Calculate distance from cluster to vertex fitted inner/outer centroids
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static float GetDistanceSquaredToCurrentVertex(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Calculate transverse impact parameter of a cluster
     * 
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static float GetImpactParameterToCurrentVertex(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Calculate the impact parameter of a cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  longitudinal the longitudinal impact parameter
     *  @param  transverse the transverse impact parameter
     */
    static void GetImpactParametersToCurrentVertex(const pandora::Cluster *const pCluster, float &longitudinal, float &transverse);

    /**
     *  @brief  Apply linear fit to inner layers of cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  innerPosition the fitted outer position
     *  @param  innerDirection the fitted outer direction
     */
    static void GetInnerVertexAndDirection(const pandora::Cluster *const pCluster, pandora::CartesianVector &innerPosition,
        pandora::CartesianVector &innerDirection);

    /**
     *  @brief  Apply linear fit to outer layers of cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  outerPosition the fitted outer position
     *  @param  outerDirection the fitted outer direction
     */
    static void GetOuterVertexAndDirection(const pandora::Cluster *const pCluster, pandora::CartesianVector &outerPosition,
        pandora::CartesianVector &outerDirection);
  
    /**
     *  @brief  Ask if cluster direction measure is correct, wrt a specified vertex
     *
     *  @param  specifiedVertex the specified vertex
     *  @param  startPosition the cluster start position
     *  @param  endPosition the cluster end position
     *  @param  direction the cluster direction measure position
     * 
     *  @return boolean
     */
    static bool IsDirectionCorrect(const pandora::CartesianVector &specifiedVertex, const pandora::CartesianVector &startPosition,
        const pandora::CartesianVector &endPosition, const pandora::CartesianVector &direction);

    /**
     *  @brief  Impact parameters
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
     *  @brief  Reset name to vertex map
     */
    static void Reset();

    /**
     *  @brief  Read the vertex helper settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    typedef std::map<std::string, pandora::CartesianVector> NameToVertexMap;

    static NameToVertexMap  m_nameToVertexMap;      ///< The name to vertex map
    static std::string      m_currentVertexName;    ///< The name of the current vertex

    static unsigned int     m_numberOfLayersToFit;  ///< The number of layers to fit for direction-finding purposes
    static float            m_longitudinalWindow;   ///< Longitudinal window used for impact parameter searches
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArVertexHelper::Reset()
{
    m_nameToVertexMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArVertexHelper::DoesVertexExist(const std::string &vertexName)
{
    return (m_nameToVertexMap.end() != m_nameToVertexMap.find(vertexName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArVertexHelper::DoesCurrentVertexExist()
{
    return DoesVertexExist(m_currentVertexName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArVertexHelper::IsForwardInZ(const pandora::Cluster *const pCluster)
{
    return LArVertexHelper::IsForwardInZ(pCluster, LArVertexHelper::GetCurrentVertex());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArVertexHelper::IsBackwardInZ(const pandora::Cluster *const pCluster)
{
    return LArVertexHelper::IsBackwardInZ(pCluster, LArVertexHelper::GetCurrentVertex());
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_HELPER_H
