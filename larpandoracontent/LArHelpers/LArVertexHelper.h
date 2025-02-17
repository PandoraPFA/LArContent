/**
 *  @file   larpandoracontent/LArHelpers/LArVertexHelper.h
 *
 *  @brief  Header file for the vertex helper class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_HELPER_H
#define LAR_VERTEX_HELPER_H 1

#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"
#include "Objects/Vertex.h"
#include "Plugins/LArTransformationPlugin.h"

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

    /**
     *  @brief  Determine if a vertex is within a detector's fiducial volume. This throws a STATUS_CODE_INVALID_PARAMETER exception if the detector
     *          is not recognised.
     *
     *  @param  pandora The Pandora instance
     *  @param  vertex The vertex to check
     *  @param  detector The string describing the detector of interest
     *                      DUNEFD HD: dune_fd_hd
     *
     *  @return true if in fiducial volume, false if not
     */
    static bool IsInFiducialVolume(const pandora::Pandora &pandora, const pandora::CartesianVector &vertex, const std::string &detector);

    /**
     *  @brief  Retrieve the true neutrino vertex position.
     *  @param  vertex The cartesian vector containing the 3D vertex position
     *  @param  x The output x coordinate
     *  @param  u The output y coordinate
     *  @param  v The output z coordinate
     */
    static void GetTrueVertexPosition(const pandora::CartesianVector &vertex, float &x, float &y, float &z);

    /**
     *  @brief  Retrieve the true neutrino vertex position.
     *
     *  @param  vertex The cartesian vector containing the 3D vertex position
     *  @param  x The output drift coordinate
     *  @param  u The output channel coordinate in the U plane
     *  @param  v The output channel coordinate in the V plane
     *  @param  w The output channel coordinate in the W plane
     */
    static void GetTrueVertexPosition(const pandora::CartesianVector &vertex, const pandora::LArTransformationPlugin *const pTransform,
        float &x, float &u, float &v, float &w);
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_HELPER_H
