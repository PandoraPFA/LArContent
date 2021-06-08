/**
 *  @file   larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h
 *
 *  @brief  Header file for the global asymmetry feature tool class.
 *
 *  $Log: $
 */
#ifndef LAR_GLOBAL_ASYMMETRY_FEATURE_TOOL_H
#define LAR_GLOBAL_ASYMMETRY_FEATURE_TOOL_H 1

#include "larpandoracontent/LArVertex/AsymmetryFeatureBaseTool.h"

namespace lar_content
{

/**
 *  @brief  GlobalAsymmetryFeatureTool class
 */
class GlobalAsymmetryFeatureTool : public AsymmetryFeatureBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    GlobalAsymmetryFeatureTool();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

private:
    /**
     *  @brief  Get the global asymmetry feature for a given view
     *
     *  @param  vertexPosition2D the vertex position projected into this view
     *  @param  slidingFitDataList the list of sliding fit data objects for this view
     *
     *  @return the global asymmetry feature
     */
    float GetAsymmetryForView(const pandora::CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList,
        const VertexSelectionBaseAlgorithm::ShowerClusterList &) const override;
};

} // namespace lar_content

#endif // #ifndef LAR_GLOBAL_ASYMMETRY_FEATURE_TOOL_H
