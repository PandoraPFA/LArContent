/**
 *  @file   larpandoracontent/LArHelpers/LArGeometryHelper.h
 *
 *  @brief  Header file for the geometry helper class.
 *
 *  $Log: $
 */
#ifndef LAR_GEOMETRY_HELPER_H
#define LAR_GEOMETRY_HELPER_H 1

#include "Pandora/PandoraEnumeratedTypes.h"
#include "Pandora/StatusCodes.h"

#include <unordered_map>

namespace pandora {class CartesianVector; class Pandora;}

namespace lar_content
{

class TwoDSlidingFitResult;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArGeometryHelper class
 */
class LArGeometryHelper
{
public:
    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     */
    static float MergeTwoPositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const float position1, const float position2);

    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  direction1 the direction in the first view
     *  @param  direction2 the direction in the second view
     */
    static pandora::CartesianVector MergeTwoDirections(const pandora::Pandora &pandora, const pandora::HitType view1,
        const pandora::HitType view2, const pandora::CartesianVector &direction1, const pandora::CartesianVector &direction2);

    /**
     *  @brief  Merge 2D positions from two views to give 2D position in third view
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 output position in the third view
     *  @param  chi-squared
     */
    static void MergeTwoPositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &position3,
        float &chiSquared);

    /**
     *  @brief  Merge 2D positions from two views to give 2D position in third view
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeTwoPositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &outputU,
        pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 2D positions for each view
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  view3 the third view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 the position in the third view
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeThreePositions(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::HitType view3, const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        const pandora::CartesianVector &position3, pandora::CartesianVector &outputU, pandora::CartesianVector &outputV,
        pandora::CartesianVector &outputW, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 2D positions for each view
     *
     *  @param  pandora the associated pandora instance
     *  @param  positionU input position in the U view
     *  @param  positionV input position in the V view
     *  @param  positionW input position in the W view
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeThreePositions(const pandora::Pandora &pandora, const pandora::CartesianVector &positionU,
        const pandora::CartesianVector &positionV, const pandora::CartesianVector &positionW, pandora::CartesianVector &outputU,
        pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, float &chiSquared);

    /**
     *  @brief  Merge 2D positions from two views to give unified 3D position
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3D output position in 3D
     *  @param  chi-squared
     */
    static void MergeTwoPositions3D(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &position3D,
        float &chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 3D position
     *
     *  @param  pandora the associated pandora instance
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  view3 the third view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 the position in the third view
     *  @param  position3D output position in 3D
     *  @param  chi-squared
     */
    static void MergeThreePositions3D(const pandora::Pandora &pandora, const pandora::HitType view1, const pandora::HitType view2,
        const pandora::HitType view3, const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        const pandora::CartesianVector &position3, pandora::CartesianVector &position3D, float &chiSquared);

    /**
     *  @brief  Project 3D position into a given 2D view
     *
     *  @param  pandora the associated pandora instance
     *  @param  position3D the position in 3D
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector ProjectPosition(const pandora::Pandora &pandora, const pandora::CartesianVector &position3D,
        const pandora::HitType view);

    /**
     *  @brief  Project 3D direction into a given 2D view
     *
     *  @param  pandora the associated pandora instance
     *  @param  direction3D the direction in 3D
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector ProjectDirection(const pandora::Pandora &pandora, const pandora::CartesianVector &direction3D,
        const pandora::HitType view);

    /**
     *  @brief  Return the wire pitch
     *
     *  @param  pandora the associated pandora instance
     *  @param  maxWirePitchDiscrepancy maximum allowed discrepancy between lar tpc wire pitch value in the w view
     */
    static float GetWireZPitch(const pandora::Pandora &pandora, const float maxWirePitchDiscrepancy = 0.01);

    /**
     *  @brief  Return the wire pitch
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     *  @param  maxWirePitchDiscrepancy maximum allowed discrepancy between lar tpc wire pitch values
     */
    static float GetWirePitch(const pandora::Pandora &pandora, const pandora::HitType view, const float maxWirePitchDiscrepancy = 0.01);

    /**
     *  @brief  Return the wire axis (vector perpendicular to the wire direction and drift direction)
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector GetWireAxis(const pandora::Pandora &pandora, const pandora::HitType view);

    /**
     *  @brief  Return the set of common daughter volumes between two 2D clusters
     *
     *  @param  Cluster 1
     *  @param  Cluster 2
     */
    static std::set<unsigned int> GetCommonDaughterVolumes (const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);
 
    /**
     *  @brief  Whether a 2D test point lies in a registered gap with the associated hit type
     *
     *  @param  pandora the associated pandora instance
     *  @param  testPoint the test point
     *  @param  hitType the hit type
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsInGap(const pandora::Pandora &pandora, const pandora::CartesianVector &testPoint2D, const pandora::HitType hitType,
        const float gapTolerance = 0.f);

    /**
     *  @brief  Whether a 3D test point lies in a registered gap with the associated hit type
     *
     *  @param  pandora the associated pandora instance
     *  @param  testPoint the test point
     *  @param  hitType the hit type
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsInGap3D(const pandora::Pandora &pandora, const pandora::CartesianVector &testPoint3D, const pandora::HitType hitType,
        const float gapTolerance = 0.f);

    /**
     *  @brief  Whether there is a gap in a cluster (described via its sliding fit result) at a specified x sampling position
     *
     *  @param  pandora the associated pandora instance
     *  @param  xSample the x sampling position
     *  @param  slidingFitResult the sliding fit result for a cluster
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsXSamplingPointInGap(const pandora::Pandora &pandora, const float xSample, const TwoDSlidingFitResult &slidingFitResult,
        const float gapTolerance = 0.f);

    /**
     *  @brief  Calculate the total distance within a given 2D region that is composed of detector gaps
     *
     *  @param  pandora the associated pandora instance
     *  @param  minZ the start position in Z
     *  @param  maxZ the end position in Z
     *  @param  hitType the hit type
     */
    static float CalculateGapDeltaZ(const pandora::Pandora &pandora, const float minZ, const float maxZ, const pandora::HitType hitType);

    /**
     *  @brief  Find the sigmaUVW value for the detector geometry
     *
     *  @param  pandora the associated pandora instance
     *  @param  maxSigmaDiscrepancy maximum allowed discrepancy between lar tpc sigmaUVW values
     */
    static float GetSigmaUVW(const pandora::Pandora &pandora, const float maxSigmaDiscrepancy = 0.01);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArGeometryHelper::GetWireZPitch(const pandora::Pandora &pandora, const float maxWirePitchWDiscrepancy)
{
    return LArGeometryHelper::GetWirePitch(pandora, pandora::TPC_VIEW_W, maxWirePitchWDiscrepancy);
}

} // namespace lar_content

#endif // #ifndef LAR_GEOMETRY_HELPER_H
