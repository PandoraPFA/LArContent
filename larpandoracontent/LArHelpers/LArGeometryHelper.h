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

namespace pandora {class CartesianVector; class Pandora; }

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

class LArPseudoLayerPlugin;
class LArTransformationPlugin;
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
     */
    static float GetWireZPitch(const pandora::Pandora &pandora);

    /**
     *  @brief  Return the wire pitch
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     */
    static float GetWirePitch(const pandora::Pandora &pandora, const pandora::HitType view);

    /**
     *  @brief  Return the wire axis (vector perpendicular to the wire direction and drift direction)
     *
     *  @param  pandora the associated pandora instance
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector GetWireAxis(const pandora::Pandora &pandora, const pandora::HitType view);

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
     *  @brief  Whether a pair of 2D test points lies in the same registered gap with the associated hit type
     *
     *  @param  pandora the associated pandora instance
     *  @param  p1 the first test point
     *  @param  p2 the second test point
     *  @param  hitType the hit type
     *  @param  gapTolerance the gap tolerance
     *
     *  @return boolean
     */
    static bool IsInGap(const pandora::Pandora &pandora, const pandora::CartesianVector &p1, const pandora::CartesianVector &p2,
        const pandora::HitType hitType, const float gapTolerance = 0.f);

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
     *  @brief  Calculate displacement of a given z-position from a specified line gap
     *
     *  @param  pLineGap  the specified line gap
     *  @param  z the given position in Z
     *
     *  @return the displacement
     */
    static float CalculateGapDisplacement(const pandora::LineGap *const pLineGap, const float &z);

    /**
     *  @brief  Calculate displacement of a given position from a specified box gap
     *
     *  @param  pBoxGap  the specified box gap
     *  @param  z the given position in Z
     *
     *  @return the displacement
     */
    static float CalculateGapDisplacement(const pandora::BoxGap *const pBoxGap, const pandora::CartesianVector &positionVector);

    /**
     *  @brief  Calculate the fractional overlap between a pair of input points and the list of registered gaps
     *
     *  @param  pandora the associated pandora instance
     *  @param  p1 the first test point
     *  @param  p2 the second test point
     *  @param  hitType the hit type
     *
     *  @return The fraction of the distance between the test points that lies within a gap.
     */
    static float CalculateOverlapFraction(const pandora::Pandora &pandora, const pandora::CartesianVector &p1,
        const pandora::CartesianVector &p2, const pandora::HitType hitType);

    /**
     *  @brief  Calculate the fractional overlap between a pair of input points and a single registered line gap
     *
     *  @param  pLineGap the input line gap
     *  @param  z1 the Z position of the first test point
     *  @param  z2 the Z position of the second test point
     *
     *  @return The fraction of the distance between the test points that lies within this line gap.
     */
    static float CalculateOverlapFraction(const pandora::LineGap *const pLineGap, const float &z1, const float &z2);

    /**
     *  @brief  Calculate the fractional overlap between a pair of input points and a single registered box gap
     *
     *  @param  pBoxGap the input box gap
     *  @param  p1 the position of the first test point
     *  @param  p2 the position of the second test point
     *
     *  @return The fraction of the distance between the test points that lies within this box gap.
     */
    static float CalculateOverlapFraction(const pandora::BoxGap *const pBoxGap, const pandora::CartesianVector &p1,
        const pandora::CartesianVector &p2);

    /**
     *  @brief  Calculate intersection between a line and a *finite* plane
     *
     *  @param  v0 the corner point on the plane
     *  @param  s1 the extent of the plane along the first edge
     *  @param  s2 the extent of the plane along the second edge
     *  @param  r0 a point on the line
     *  @param  traj the direction vector of the line
     *
     *  @return The intersection point if it exists
     */
    static pandora::CartesianVector CalculateIntersection(const pandora::CartesianVector &v0, const pandora::CartesianVector &s1,
        const pandora::CartesianVector &s2, const pandora::CartesianVector &r0, const pandora::CartesianVector &traj);

    /**
     *  @brief  Get the LArPseudoLayerPlugin registered with a specified pandora instance
     *
     *  @param  pandora the associated pandora instance
     *
     *  @return the address of the relevant LArPseudoLayerPlugin
     */
    static const LArPseudoLayerPlugin *GetLArPseudoLayerPlugin(const pandora::Pandora &pandora);

    /**
     *  @brief  Get the LArTransformationPlugin registered with a specified pandora instance
     *
     *  @param  pandora the associated pandora instance
     *
     *  @return the address of the relevant LArTransformationPlugin
     */
    static const LArTransformationPlugin *GetLArTransformationPlugin(const pandora::Pandora &pandora);

    /**
     *  @brief  Set the LArPseudoLayerPlugin for a given pandora instance
     *
     *  @param  pandora the pandora instance
     *  @param  pLArPseudoLayerPlugin the address of the LArPseudoLayerPlugin
     */
    static pandora::StatusCode SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, const LArPseudoLayerPlugin *const pLArPseudoLayerPlugin);

    /**
     *  @brief  Set the LArTransformationPlugin for a given pandora instance
     *
     *  @param  pandora the pandora instance
     *  @param  pLArTransformationPlugin the address of the LArTransformationPlugin
     */
    static pandora::StatusCode SetLArTransformationPlugin(const pandora::Pandora &pandora, const LArTransformationPlugin *const pLArTransformationPlugin);

private:
    typedef std::unordered_map<const pandora::Pandora*, const LArPseudoLayerPlugin*>    PseudoLayerInstanceMap;
    typedef std::unordered_map<const pandora::Pandora*, const LArTransformationPlugin*> TransformationInstanceMap;
    static PseudoLayerInstanceMap       m_pseudolayerInstanceMap;       ///< The pseudolayer instance map
    static TransformationInstanceMap    m_transformationInstanceMap;    ///< The transformation instance map
};

} // namespace lar_content

#endif // #ifndef LAR_GEOMETRY_HELPER_H
