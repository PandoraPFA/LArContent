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
     *  @brief  Merge 2D positions from two views to give 2D position in third view
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *
     *  @param  position3 output position in the U view
     *  @param  chi-squared
     */
    static void MergeTwoPositions(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        pandora::CartesianVector &position3, float& chiSquared);

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
    static void MergeTwoPositions(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2,
        pandora::CartesianVector &outputU, pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, 
        float& chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 2D positions for each view
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
    static void MergeThreePositions(const pandora::HitType view1, const pandora::HitType view2, const pandora::HitType view3, 
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, const pandora::CartesianVector &position3,  
        pandora::CartesianVector &outputU, pandora::CartesianVector &outputV, pandora::CartesianVector &outputW, 
        float& chiSquared);

    /**
     *  @brief  Merge 2D positions from three views to give unified 2D positions for each view
     *
     *  @param  positionU input position in the U view
     *  @param  positionV input position in the V view
     *  @param  positionW input position in the W view 
     *  @param  positionU output position in the U view
     *  @param  positionV output position in the V view
     *  @param  positionW output position in the W view
     *  @param  chi-squared
     */
    static void MergeThreePositions(const pandora::CartesianVector &positionU, const pandora::CartesianVector &positionV, 
        const pandora::CartesianVector &positionW, pandora::CartesianVector &outputU, pandora::CartesianVector &outputV, 
        pandora::CartesianVector &outputW, float& chiSquared);


    /**
     *  @brief  Merge 2D positions from two views to give unified 3D position
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *
     *  @param  position3D output position in 3D
     */
    static void MergeTwoPositions3D(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &position3D);

    /**
     *  @brief  Merge 2D positions from three views to give unified 3D position
     *
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  view3 the third view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     *  @param  position3 the position in the third view
     *
     *  @param  position3D output position in 3D
     */
    static void MergeThreePositions3D(const pandora::HitType view1, const pandora::HitType view2, const pandora::HitType view3,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, const pandora::CartesianVector &position3, 
        pandora::CartesianVector &position3D);


    /**
     *  @brief  Project 3D position into a given 2D view
     *
     *  @param  position3D the position in 3D
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector ProjectPosition(const pandora::CartesianVector &position3D, const pandora::HitType view);

    /**
     *  @brief  Project 3D direction into a given 2D view
     *
     *  @param  direction3D the direction in 3D
     *  @param  view the 2D projection
     */
    static pandora::CartesianVector ProjectDirection(const pandora::CartesianVector &direction3D, const pandora::HitType view);


    /** 
     *  @brief  Transform from (U,V) to W position 
     *
     *  @param U the U position
     *  @param V the V position  
     */
     static float UVtoW(const float& U, const float& V);

    /** 
     *  @brief  Transform from (V,W) to U position 
     *
     *  @param V the V position
     *  @param W the W position  
     */
     static float VWtoU(const float& V, const float& W);

    /** 
     *  @brief  Transform from (W,U) to V position 
     *
     *  @param W the W position
     *  @param U the U position  
     */
     static float WUtoV(const float& W, const float& U); 

    /** 
     *  @brief  Transform from (V,U) to W position 
     * 
     *  @param V the V position  
     *  @param U the U position
     */
     static float VUtoW(const float& V, const float& U);

    /** 
     *  @brief  Transform from (W,V) to U position 
     *
     *  @param W the W position  
     *  @param V the V position
     */
     static float WVtoU(const float& W, const float& V);

    /** 
     *  @brief  Transform from (U,W) to V position 
     *
     *  @param U the U position  
     *  @param W the W position
     */
     static float UWtoV(const float& U, const float& W); 


    /** 
     *  @brief  Transform from (U,V) to Y position 
     * 
     *              Y * sin(thetaU+thetaV) = U * cos(thetaV) + V * cos(thetaU) - H/2 * sin(thetaU-thetaV)
     *
     *  @param U the U position
     *  @param V the V position  
     */
    static float UVtoY(const float& U, const float& V);

    /** 
     *  @brief  Transform from (U,V) to Z position 
     * 
     *              Z * sin(thetaU+thetaV) = U * sin(thetaV) + V * sin(thetaU) - H * sin(thetaU) * sin(thetaV)
     *
     *  @param U the U position
     *  @param V the V position  
     */
     static float UVtoZ(const float& U, const float& V);  

    /** 
     *  @brief  Transform from (Y,Z) to U position
     *
     *              U = Z * cos(thetaU) + ( Y + H/2 ) * sin(thetaU)
     *   
     *  @param Y the Y position   
     *  @param Z the Z position   
     */
     static float YZtoU(const float& Y, const float& Z);

    /** 
     *  @brief  Transform from (Y,Z) to V position
     *
     *              V = Z * cos(thetaV) - ( Y - H/2 ) * sin(thetaV)
     *   
     *  @param Y the Y position   
     *  @param Z the Z position   
     */
     static float YZtoV(const float& Y, const float& Z);

    /** 
     *  @brief  Transform from (pU,pV) to pW direction 
     *
     *  @param pU the pU direction
     *  @param pV the pV direction  
     */
     static float PUPVtoPW(const float& pU, const float& pV);

    /** 
     *  @brief  Transform from (pV,pW) to pU direction 
     *
     *  @param pV the pV direction
     *  @param pW the pW direction  
     */
     static float PVPWtoPU(const float& pV, const float& pW);

    /** 
     *  @brief  Transform from (pW,pU) to pV direction 
     *
     *  @param pW the pW direction
     *  @param pU the pU direction  
     */
     static float PWPUtoPV(const float& pW, const float& pU); 

    /** 
     *  @brief  Transform from (pV,pU) to pW direction 
     * 
     *  @param pV the pV direction  
     *  @param pU the pU direction
     */
     static float PVPUtoPW(const float& pV, const float& pU);

    /** 
     *  @brief  Transform from (pW,pV) to pU direction 
     *
     *  @param pW the pW direction  
     *  @param pV the pV direction
     */
     static float PWPVtoPU(const float& pW, const float& pV);

    /** 
     *  @brief  Transform from (pU,pW) to pV direction 
     *
     *  @param pU the pU direction  
     *  @param pW the pW direction
     */
     static float PUPWtoPV(const float& pU, const float& pW); 

    /** 
     *  @brief  Transform from (pY,pZ) to pU direction
     *
     *              pU = pZ * cos(thetaU) + pY * sin(thetaU)  
     *   
     *  @param pU the U component   
     *  @param pV the V component   
     */
    static float PYPZtoPU(const float& pY, const float& pZ); 

    /** 
     *  @brief  Transform from (pY,pZ) to pV direction
     *
     *              pV = pZ * cos(thetaV) - pY * sin(thetaV)
     *   
     *  @param pU the U component   
     *  @param pV the V component   
     */
    static float PYPZtoPV(const float& pY, const float& pZ);


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
