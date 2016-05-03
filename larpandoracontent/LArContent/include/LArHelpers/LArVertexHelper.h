/**
 *  @file   LArContent/include/LArHelpers/LArVertexHelper.h
 *
 *  @brief  Header file for the vertex helper class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_HELPER_H
#define LAR_VERTEX_HELPER_H 1

#include "Objects/Cluster.h"
#include "Objects/Vertex.h"

namespace lar_content
{

/**
 *  @brief  LArVertexHelper class
 */
class LArVertexHelper
{
public:
    /**
     *  ClusterDirection enumeration
     */
    enum ClusterDirection
    {
        DIRECTION_FORWARD_IN_Z,
        DIRECTION_BACKWARD_IN_Z,
        DIRECTION_UNKNOWN
    };

    /**
     *  @brief  Get the direction of the cluster in z, using a projection of the provided vertex
     * 
     *  @param  pandora the pandora instance
     *  @param  pVertex the address of the vertex
     *  @param  pCluster the address of the cluster
     *  @param  tanAngle look for vertex inside triangle with apex shifted along the cluster length
     *  @param  apexShift look for vertex inside triangle with apex shifted along the cluster length
     * 
     *  @return the cluster direction in z
     */
    static ClusterDirection GetClusterDirectionInZ(const pandora::Pandora &pandora, const pandora::Vertex *const pVertex,
        const pandora::Cluster *const pCluster, const float tanAngle, const float apexShift);
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_HELPER_H
