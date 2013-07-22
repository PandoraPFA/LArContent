/**
 *  @file   LArContent/include/LArHelpers/LArGeometryHelper.h
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
     *  @brief  Merge two views (U,V) to give a third view (Z). 
     *          U and V views are inclined at angles thetaU and thetaV (both positive)
     *          and H is the height of the detector
     * 
     *          Transformation from (Y,Z) to (U,V) coordinates:
     *            U = Z * cos(thetaU) + ( Y + H/2 ) * sin(thetaU)
     *            V = Z * cos(thetaV) - ( Y - H/2 ) * sin(thetaV)
     *
     *          (Eliminating Y...)
     *         => Z * sin(thetaU+thetaV) = U * sin(thetaV) + V * sin(thetaU) - H * sin(thetaU) * sin(thetaV)
     *                                                       
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     */
    static float MergeTwoPositions(const pandora::HitType view1, const pandora::HitType view2, const float position1, const float position2);

    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).  
     *          U and V views are inclined at angles thetaU and thetaV (both positive)
     * 
     *          Transformation from (pY,pZ) to (pU,pV) coordinates:
     *            pU = pZ * cos(thetaU) + pY * sin(thetaU)
     *            pV = pZ * cos(thetaV) - pY * sin(thetaV)
     *
     *          (Eliminating pY...)
     *         => pZ * sin(thetaU+thetaV) = pU * sin(thetaV) + pV * sin(thetaU)            
     * 
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  direction1 the direction in the first view
     *  @param  direction2 the direction in the second view
     */
    static pandora::CartesianVector MergeTwoDirections(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &direction1, const pandora::CartesianVector &direction2);

    /**
     *  @brief  Merge two views (1,2)
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *
     *  @param  position3 output position in the U view
     *  @param  chi-squared
     */
    static void MergeTwoViews(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        pandora::CartesianVector &position3, float& chiSquared);

    /**
     *  @brief  Merge two views (1,2)
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
    static void MergeTwoViews(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        pandora::CartesianVector &outputU, pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, 
        float& chiSquared);

    /**
     *  @brief  Merge three views (1,2,3)
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  view3 the third view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 the position in the third view
     *
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view 
     *  @param  chi-squared
     */
    static void MergeThreeViews(const pandora::HitType view1, const pandora::HitType view2, const pandora::HitType view3, 
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, const pandora::CartesianVector &position3,  
        pandora::CartesianVector &outputU, pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, 
        float& chiSquared);

    /**
     *  @brief  Merge three views (U,V,W). 
     *
     *  @param  positionU input position in the U view
     *  @param  positionV input position in the V view
     *  @param  positionW input position in the W view 
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeThreeViews(const pandora::CartesianVector &positionU, const pandora::CartesianVector &positionV, 
        const pandora::CartesianVector &positionW, pandora::CartesianVector &outputU, pandora::CartesianVector &outputV, 
        pandora::CartesianVector &outputW, float& chiSquared);

    /** 
     *  @brief  Transform from (U,V) to W coordinate 
     *
     *  @param U the U coordinate
     *  @param V the V coordinate  
     */

     static float UVtoW(const float& U, const float& V);

    /** 
     *  @brief  Transform from (V,W) to U coordinate 
     *
     *  @param V the V coordinate
     *  @param W the W coordinate  
     */

     static float VWtoU(const float& V, const float& W);

    /** 
     *  @brief  Transform from (W,U) to V coordinate 
     *
     *  @param W the W coordinate
     *  @param U the U coordinate  
     */

     static float WUtoV(const float& W, const float& U); 

    /** 
     *  @brief  Transform from (V,U) to W coordinate 
     * 
     *  @param V the V coordinate  
     *  @param U the U coordinate
     */

     static float VUtoW(const float& V, const float& U);

    /** 
     *  @brief  Transform from (W,V) to U coordinate 
     *
     *  @param W the W coordinate  
     *  @param V the V coordinate
     */

     static float WVtoU(const float& W, const float& V);

    /** 
     *  @brief  Transform from (U,W) to V coordinate 
     *
     *  @param U the U coordinate  
     *  @param W the W coordinate
     */

     static float UWtoV(const float& U, const float& W); 


    /** 
     *  @brief  Transform from (U,V) to Y coordinate 
     * 
     *              Y * sin(thetaU+thetaV) = U * cos(thetaV) + V * cos(thetaU) - H/2 * sin(thetaU-thetaV)
     *
     *  @param U the U coordinate
     *  @param V the V coordinate  
     */

    static float UVtoY(const float& U, const float& V);

    /** 
     *  @brief  Transform from (U,V) to Z coordinate 
     * 
     *              Z * sin(thetaU+thetaV) = U * sin(thetaV) + V * sin(thetaU) - H * sin(thetaU) * sin(thetaV)
     *
     *  @param U the U coordinate
     *  @param V the V coordinate  
     */

     static float UVtoZ(const float& U, const float& V);  

    /** 
     *  @brief  Transform from (Y,Z) to U coordinate
     *
     *              U = Z * cos(thetaU) + ( Y + H/2 ) * sin(thetaU)
     *   
     *  @param Y the Y coordinate   
     *  @param Z the Z coordinate   
     */

     static float YZtoU(const float& Y, const float& Z);

    /** 
     *  @brief  Transform from (Y,Z) to V coordinate
     *
     *              V = Z * cos(thetaV) - ( Y - H/2 ) * sin(thetaV)
     *   
     *  @param Y the Y coordinate   
     *  @param Z the Z coordinate   
     */

    static float YZtoV(const float& Y, const float& Z);

    /**
     *  @brief  Read the vertex helper settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */

    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    static float    m_thetaU;            // inclination of U wires (radians)
    static float    m_thetaV;            // inclination of V wires (radians)

    static float    m_sinUminusV;        // sin(thetaU-thetaV)
    static float    m_sinUplusV;         // sin(thetaU+thetaV)
    static float    m_sinU;              // sin(thetaU)
    static float    m_sinV;              // sin(thetaV)
    static float    m_cosU;              // cos(thetaU)
    static float    m_cosV;              // cos(thetaV)
    static float    m_H;                 // height (cm)

    static float    m_sigmaUVW;          // resolution (cm), for calculation of chi2
};

} // namespace lar

#endif // #ifndef LAR_GEOMETRY_HELPER_H
