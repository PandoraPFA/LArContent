/**
 *  @file   larpandoracontent/LArHelpers/LArStitchingHelper.h
 *
 *  @brief  Header file for the helper class for multiple drift volumes.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_HELPER_H
#define LAR_STITCHING_HELPER_H 1

#include "Objects/ParticleFlowObject.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

namespace lar_content
{

/**
 *  @brief  LArStitchingHelper class
 */
class LArStitchingHelper
{
public:
   /**
     *  @brief Find neighbouring drift volume to a specified drift volume
     *
     *  @param pandora  the pandora stitching instance
     *  @param inputVolume  the specified drift volume
     *  @param checkPositive  look in higher (lower) x positions if this is set to true (false)
     *
     *  @return volume information block
     */
    static const VolumeInfo &FindClosestVolume(const pandora::Pandora &pandora, const VolumeInfo &inputVolume, const bool checkPositive);

    /**
     *  @brief Check that a pair of drift volumes are adjacent to each other
     *
     *  @param firstVolume the first drift volume
     *  @param secondVolume the second drift volume
     *
     *  @return boolean
     */
    static bool CanVolumesBeStitched(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume);

    /**
     *  @brief Check that a pair of drift volumes are adjacent to each other
     *
     *  @param firstVolume the first drift volume
     *  @param secondVolume the second drift volume
     *
     *  @return boolean
     */
    static bool AreVolumesAdjacent(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume);

    /**
     *  @brief Check that a pair of drift volumes are adjacent to each other
     *
     *  @param pandora  the pandora stitching instance
     *  @param firstVolume  the first drift volume
     *  @param secondVolume  the second drift volume
     *
     *  @return boolean
     */
    static bool AreVolumesAdjacent(const pandora::Pandora &pandora, const VolumeInfo &firstVolume, const VolumeInfo &secondVolume);

    /**
     *  @brief Determine centre in X at the boundary between a pair of drift volumes
     *
     *  @param firstVolume the first drift volume
     *  @param secondVolume the second drift volume
     *
     *  @return Boundary X centre
     */
    static float GetVolumeBoundaryCenterX(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume);

    /**
     *  @brief Determine width in X at the boundary between a pair of drift volumes
     *
     *  @param firstVolume the first drift volume
     *  @param secondVolume the second drift volume
     *
     *  @return Boundary X width
     */
    static float GetVolumeBoundaryWidthX(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume);

    /**
     *  @brief Calculate distance between central positions of a pair of drift volumes
     *
     *  @param firstVolume the first drift volume
     *  @param secondVolume the second drift volume
     *
     *  @return the distance as a float
     */
    static float GetVolumeDisplacement(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume);

    /**
     *  @brief  Given a pair of pointing clusters, find the pair of vertices with smallest yz-separation
     *
     *  @param  driftVolume1 the first drift volume
     *  @param  driftVolume2 the second drift volume
     *  @param  pointingCluster1 the pointing cluster in the first drift volume
     *  @param  pointingCluster2 the pointing cluster in the second drift volume
     *  @param  closestVertex1 to receive the relevant vertex from the first pointing cluster
     *  @param  closestVertex2 to receive the relevant vertex from the second pointing cluster
     */
    static void GetClosestVertices(const VolumeInfo &driftVolume1, const VolumeInfo &driftVolume2,
        const LArPointingCluster &pointingCluster1, const LArPointingCluster &pointingCluster2,
        LArPointingCluster::Vertex &closestVertex1, LArPointingCluster::Vertex &closestVertex2);

    /**
     *  @brief  Calculate X0 for a pair of vertices
     *
     *  @param  firstVolume the first drift volume
     *  @param  secondVolume the second drift volume
     *  @param  firstVertex the relevant vertex from the first pointing cluster
     *  @param  secondVertex the relevant vertex from the second pointing cluster
     *
     *  @return X0 value for this pair of vertices
     */
    static float CalculateX0(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume,
        const LArPointingCluster::Vertex &firstVertex, const LArPointingCluster::Vertex &secondVertex);

    /**
     *  @brief  Apply the X0 correction to an input position vector
     *
     *  @param  driftVolume  the drift volume
     *  @param  x0  the x0 value
     *  @param  inputPosition  the input uncorrected position vector
     *
     *  @return the output corrected position vector
     */
    static pandora::CartesianVector GetCorrectedPosition(const VolumeInfo &driftVolume, const float x0,
        const pandora::CartesianVector &inputPosition);

};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_HELPER_H
