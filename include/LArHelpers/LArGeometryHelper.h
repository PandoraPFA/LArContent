/**
 *  @file   LArContent/include/LArHelpers/LArGeometryHelper.h
 * 
 *  @brief  Header file for the geometry helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_GEOMETRY_HELPER_H
#define LAR_GEOMETRY_HELPER_H 1

#include "Api/PandoraApi.h"

#include "Objects/CartesianVector.h"

namespace lar
{

class LArPseudoLayerCalculator;
class LArTransformationCalculator;

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
     *  @param  view1 the first view
     *  @param  view2 the second view
     *  @param  position1 the position in the first view
     *  @param  position2 the position in the second view
     */
    static float MergeTwoPositions(const pandora::HitType view1, const pandora::HitType view2, const float position1, const float position2);

    /**
     *  @brief  Merge two views (U,V) to give a third view (Z).
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
     *  @param  position3 output position in the third view
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
     *
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
     *  @param  chi-squared
     */
    static void MergeTwoPositions3D(const pandora::HitType view1, const pandora::HitType view2,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, 
        pandora::CartesianVector &position3D, float &chiSquared);

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
     *  @param  chi-squared
     */
    static void MergeThreePositions3D(const pandora::HitType view1, const pandora::HitType view2, const pandora::HitType view3,
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, const pandora::CartesianVector &position3, 
        pandora::CartesianVector &position3D, float &chiSquared);

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
     *  @brief  Set the registered LArPseudoLayerCalculator (derived version of registered geometry helper PseudoLayerCalculator)
     *
     *  @param  pLArPseudoLayerCalculator the address of the LArPseudoLayerCalculator
     */
    static pandora::StatusCode SetLArPseudoLayerCalculator(const LArPseudoLayerCalculator *pLArPseudoLayerCalculator);

    /**
     *  @brief  Get the registered LArPseudoLayerCalculator (derived version of registered geometry helper PseudoLayerCalculator)
     *
     *  @return the address of the LArPseudoLayerCalculator
     */
    static const LArPseudoLayerCalculator *GetLArPseudoLayerCalculator();

    /**
     *  @brief  Set the registered LArTransformationCalculator
     *
     *  @param  pLArTransformationCalculator the address of the LArTransformationCalculator
     */
    static pandora::StatusCode SetLArTransformationCalculator(const LArTransformationCalculator *pLArTransformationCalculator);

    /**
     *  @brief  Get the registered LArTransformationCalculator
     *
     *  @return the address of the LArTransformationCalculator
     */
    static const LArTransformationCalculator *GetLArTransformationCalculator();

    /**
     *  @brief  Read the vertex helper settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    static const LArPseudoLayerCalculator     *m_pLArPseudoLayerCalculator;       ///< Address of the lar pseudolayer calculator
    static const LArTransformationCalculator  *m_pLArTransformationCalculator;    ///< Address of the lar transformation calculator
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArGeometryHelper::SetLArPseudoLayerCalculator(const LArPseudoLayerCalculator *pLArPseudoLayerCalculator)
{
    if (NULL != m_pLArPseudoLayerCalculator)
        return pandora::STATUS_CODE_ALREADY_INITIALIZED;

    m_pLArPseudoLayerCalculator = pLArPseudoLayerCalculator;
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArPseudoLayerCalculator *LArGeometryHelper::GetLArPseudoLayerCalculator()
{
    if (NULL == m_pLArPseudoLayerCalculator)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_pLArPseudoLayerCalculator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArGeometryHelper::SetLArTransformationCalculator(const LArTransformationCalculator *pLArTransformationCalculator)
{
    if (NULL != m_pLArTransformationCalculator)
        return pandora::STATUS_CODE_ALREADY_INITIALIZED;

    m_pLArTransformationCalculator = pLArTransformationCalculator;
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArTransformationCalculator *LArGeometryHelper::GetLArTransformationCalculator()
{
    if (NULL == m_pLArTransformationCalculator)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_pLArTransformationCalculator;
}

} // namespace lar

#endif // #ifndef LAR_GEOMETRY_HELPER_H
