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

namespace lar_content
{

class LArPseudoLayerPlugin;
class LArTransformationPlugin;

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
        const pandora::CartesianVector &position1, const pandora::CartesianVector &position2, pandora::CartesianVector &position3D, float &chiSquared);

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
     *
     *  @return the wire pitch 
     */
    static float GetWireZPitch(const pandora::Pandora &pandora);

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
    static pandora::StatusCode SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, const LArPseudoLayerPlugin *pLArPseudoLayerPlugin);
  
    /**
     *  @brief  Set the LArTransformationPlugin for a given pandora instance
     *
     *  @param  pandora the pandora instance
     *  @param  pLArTransformationPlugin the address of the LArTransformationPlugin
     */
    static pandora::StatusCode SetLArTransformationPlugin(const pandora::Pandora &pandora, const LArTransformationPlugin *pLArTransformationPlugin);

private:
    typedef std::map<const pandora::Pandora*, const LArPseudoLayerPlugin*>    PseudoLayerInstanceMap;
    typedef std::map<const pandora::Pandora*, const LArTransformationPlugin*> TransformationInstanceMap;
    static PseudoLayerInstanceMap       m_pseudolayerInstanceMap;       ///< The pseudolayer instance map 
    static TransformationInstanceMap    m_transformationInstanceMap;    ///< The transformation instance map
};

} // namespace lar_content

#endif // #ifndef LAR_GEOMETRY_HELPER_H
