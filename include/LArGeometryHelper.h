/**
 *  @file   LArGeometryHelper.h
 * 
 *  @brief  Header file for the geometry helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_GEOMETRY_HELPER_H
#define LAR_GEOMETRY_HELPER_H 1

#include "Objects/CartesianVector.h"

namespace lar
{

/**
 *  @brief  LArGeometryHelper class
 */
class LArGeometryHelper
{
public:
    /**
     *  @brief  Merge two views (U,V) to give a third view (Z). U and V views are inclined at angles thetaU and thetaV (both positive)
     *          and H is the height of the detector
     * 
     *          (1) Tranformation from (Y,Z) to (U,V) coordinates:
     *              U = Z * cos(thetaU) + ( Y + H/2 ) * sin(thetaU)
     *              V = Z * cos(thetaV) - ( Y - H/2 ) * sin(thetaV)
     * 
     *          (2) Transformation between (U,V,Z) coordinates:
     *          (a) U,V -> Z:           U * sin(thetaV) + V * sin(thetaU) - H * sin(thetaU) * sin(thetaV)
     *                              Z = -----------------------------------------------------------------
     *                                                   sin(thetaU+thetaV)
     * 
     *          (b) V,Z -> U:          Z * sin(thetaU+thetaV) - V * sin(thetaU) + H * sin(thetaU) * sin(thetaV) 
     *                              U = ------------------------------------------------------------------------
     *                                                       sin(thetaV)
     * 
     *          (c) U,Z -> V:         Z * sin(thetaU+thetaV) - U * sin(thetaV) + H * sin(thetaU) * sin(thetaV) 
     *                              V = ------------------------------------------------------------------------
     *                                                       sin(thetaU)
     * 
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     */
    static float MergeTwoPositions(const pandora::HitType view1, const pandora::HitType view2, const float position1, const float position2);

    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).  U and V views are inclined at angles thetaU and thetaV (both positive)
     * 
     *          (1) Tranformation from (pY,pZ) to (pU,pV) coordinates:
     *              pU = pZ * cos(thetaU) + pY * sin(thetaU)
     *              pV = pZ * cos(thetaV) - pY * sin(thetaV)
     * 
     *          (2) Transformation between (pU,pV,pZ) coordinates:
     *          (a) pU,pV -> pZ:         pU * sin(thetaV) + pV * sin(thetaU)
     *                              pZ = -----------------------------------
     *                                           sin(thetaU+thetaV)
     * 
     *          (b) pV,pZ -> pU:         pZ * sin(thetaU+thetaV) - pV * sin(thetaU)
     *                              pU = ------------------------------------------
     *                                              sin(thetaV)
     * 
     *          (c) pU,pZ -> pV:         pZ * sin(thetaU+thetaV) - pU * sin(thetaV)
     *                              pV = ------------------------------------------
     *                                              sin(thetaU)
     * 
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  direction1 the direction in the first view
     *  @param  direction2 the direction in the second view
     */
    static pandora::CartesianVector MergeTwoDirections(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &direction1, const pandora::CartesianVector &direction2);

    /**
     *  @brief  Read the vertex helper settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    static float    m_thetaU;            // inclination of U wires (radians)
    static float    m_thetaV;            // inclination of V wires (radians)

    static float    m_sinUV;             // sin(thetaU+thetaV)
    static float    m_sinU;              // sin(thetaU)
    static float    m_sinV;              // sin(thetaV)
    static float    m_H;                 // height (cm)
};

} // namespace lar

#endif // #ifndef LAR_GEOMETRY_HELPER_H
