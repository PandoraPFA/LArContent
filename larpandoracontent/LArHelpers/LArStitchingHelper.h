/**
 *  @file   larpandoracontent/LArHelpers/LArStitchingHelper.h
 *
 *  @brief  Header file for the helper class for multiple drift volumes.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_HELPER_H
#define LAR_STITCHING_HELPER_H 1

#include "Geometry/LArTPC.h"

#include "Objects/ParticleFlowObject.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

namespace lar_content
{

/**
 *  @brief  LArStitchingHelper class
 */
class LArStitchingHelper
{
public:
   /**
     *  @brief  Find closest tpc to a specified input tpc
     *
     *  @param  pandora the pandora stitching instance
     *  @param  inputTPC the specified drift volume
     *  @param  checkPositive look in higher (lower) x positions if this is set to true (false)
     *
     *  @return the closest tpc
     */
    static const pandora::LArTPC &FindClosestTPC(const pandora::Pandora &pandora, const pandora::LArTPC &inputTPC, const bool checkPositive);

    /**
     *  @brief  Whether particles from a given pair of tpcs can be stitched together
     *
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *
     *  @return boolean
     */
    static bool CanTPCsBeStitched(const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC);

    /**
     *  @brief  Whether a pair of drift volumes are adjacent to each other
     *
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *
     *  @return boolean
     */
    static bool AreTPCsAdjacent(const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC);

    /**
     *  @brief  Whether a pair of drift volumes are adjacent to each other
     *
     *  @param  pandora  the pandora stitching instance
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *
     *  @return boolean
     */
    static bool AreTPCsAdjacent(const pandora::Pandora &pandora, const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC);

    /**
     *  @brief  Determine centre in X at the boundary between a pair of tpcs
     *
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *
     *  @return boundary X centre
     */
    static float GetTPCBoundaryCenterX(const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC);

    /**
     *  @brief  Determine width in X at the boundary between a pair of tpcs
     *
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *
     *  @return boundary X width
     */
    static float GetTPCBoundaryWidthX(const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC);

    /**
     *  @brief  Calculate distance between central positions of a pair of tpcs
     *
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *
     *  @return the distance
     */
    static float GetTPCDisplacement(const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC);

    /**
     *  @brief  Given a pair of pointing clusters, find the pair of vertices with smallest yz-separation
     *
     *  @param  larTPC1 the first tpc
     *  @param  larTPC2 the second tpc
     *  @param  pointingCluster1 the pointing cluster in the first tpc
     *  @param  pointingCluster2 the pointing cluster in the second tpc
     *  @param  closestVertex1 to receive the relevant vertex from the first pointing cluster
     *  @param  closestVertex2 to receive the relevant vertex from the second pointing cluster
     */
    static void GetClosestVertices(const pandora::LArTPC &larTPC1, const pandora::LArTPC &larTPC2,
        const LArPointingCluster &pointingCluster1, const LArPointingCluster &pointingCluster2,
        LArPointingCluster::Vertex &closestVertex1, LArPointingCluster::Vertex &closestVertex2);

    /**
     *  @brief  Calculate X0 for a pair of vertices
     *
     *  @param  firstTPC the first tpc
     *  @param  secondTPC the second tpc
     *  @param  firstVertex the relevant vertex from the first pointing cluster
     *  @param  secondVertex the relevant vertex from the second pointing cluster
     *
     *  @return X0 value for this pair of vertices
     */
    static float CalculateX0(const pandora::LArTPC &firstTPC, const pandora::LArTPC &secondTPC,
        const LArPointingCluster::Vertex &firstVertex, const LArPointingCluster::Vertex &secondVertex);

    /**
     *  @brief  Apply the X0 correction to an input position vector
     *
     *  @param  larTPC the tpc
     *  @param  x0  the x0 value
     *  @param  inputPosition  the input uncorrected position vector
     *
     *  @return the output corrected position vector
     */
    static pandora::CartesianVector GetCorrectedPosition(const pandora::LArTPC &larTPC, const float x0,
        const pandora::CartesianVector &inputPosition);

    /**
     *  @brief  Sort tpcs by central positions
     *
     *  @param  pLhs address of first tpc
     *  @param  pRhs address of second tpc
     */
    static bool SortTPCs(const pandora::LArTPC *const pLhs, const pandora::LArTPC *const pRhs);
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_HELPER_H
