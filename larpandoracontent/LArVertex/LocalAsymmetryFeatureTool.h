/**
 *  @file   larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h
 *
 *  @brief  Header file for the local asymmetry feature tool class.
 *
 *  $Log: $
 */
#ifndef LAR_LOCAL_ASYMMETRY_FEATURE_TOOL_H
#define LAR_LOCAL_ASYMMETRY_FEATURE_TOOL_H 1

#include "larpandoracontent/LArVertex/AsymmetryFeatureBaseTool.h"

namespace lar_content
{

/**
 *  @brief  LocalAsymmetryFeatureTool class
 */
class LocalAsymmetryFeatureTool : public AsymmetryFeatureBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    LocalAsymmetryFeatureTool();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    /**
     *  @brief  Get the local asymmetry feature in a given view
     *
     *  @param  vertexPosition2D the vertex position projected into this view
     *  @param  slidingFitDataList the list of sliding fit data objects in this view
     *
     *  @return the local asymmetry feature
     */
    float GetAsymmetryForView(
	const pandora::CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList,
	const VertexSelectionBaseAlgorithm::ShowerClusterList &) const override;

    /**
     *  @brief  Check whether a cluster's direction agrees with the current weighted direction
     *
     *  @param  weightedDirectionSum the current weighted direction
     *  @param  clusterDirection the cluster's direction
     *
     *  @return boolean
     */
    bool CheckAngle(const pandora::CartesianVector &weightedDirectionSum, const pandora::CartesianVector &clusterDirection) const;

    float m_minAsymmetryCosAngle;         ///< The min opening angle cosine used to determine viability of asymmetry score
    unsigned int m_maxAsymmetryNClusters; ///< The max number of associated clusters to calculate the asymmetry
};

} // namespace lar_content

#endif // #ifndef LAR_LOCAL_ASYMMETRY_FEATURE_TOOL_H
