/**
 *  @file   larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h
 *
 *  @brief  Header file for the global asymmetry feature tool class.
 *
 *  $Log: $
 */
#ifndef LAR_GLOBAL_ASYMMETRY_FEATURE_TOOL_H
#define LAR_GLOBAL_ASYMMETRY_FEATURE_TOOL_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  GlobalAsymmetryFeatureTool class
 */
class GlobalAsymmetryFeatureTool : public VertexSelectionBaseAlgorithm::VertexFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    GlobalAsymmetryFeatureTool();

    /**
     *  @brief  Run the tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pVertex address of the vertex
     *  @param  slidingFitDataListMap map of the sliding fit data lists
     *
     *  @return the global asymmetry feature
     */
    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm, const pandora::Vertex *const pVertex,
        const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap, const VertexSelectionBaseAlgorithm::ClusterListMap &,
        const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the global asymmetry feature for a given view
     *
     *  @param  vertexPosition2D the vertex position projected into this view
     *  @param  slidingFitDataList the list of sliding fit data objects for this view
     *
     *  @return the global asymmetry feature
     */
    float GetGlobalAsymmetryForView(
        const pandora::CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList) const;

    /**
     *  @brief  Increment the asymmetry parameters
     *
     *  @param  weight the weight to assign to this vector
     *  @param  clusterDirection the direction of the cluster
     *  @param  localWeightedDirectionSum the current energy-weighted local cluster direction vector
     */
    void IncrementAsymmetryParameters(
        const float weight, const pandora::CartesianVector &clusterDirection, pandora::CartesianVector &localWeightedDirectionSum) const;

    /**
     *  @brief  Calculate the global asymmetry feature
     *
     *  @param  useEnergyMetrics whether to use energy-based metrics instead of hit-counting-based metrics
     *  @param  vertexPosition2D the vertex position in this view
     *  @param  slidingFitDataList the list of sliding fit data objects
     *  @param  localWeightedDirectionSum the local event axis
     *
     *  @return the global asymmetry feature
     */
    float CalculateGlobalAsymmetry(const bool useEnergyMetrics, const pandora::CartesianVector &vertexPosition2D,
        const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList, const pandora::CartesianVector &localWeightedDirectionSum) const;

    float m_maxAsymmetryDistance; ///< The max distance between cluster (any hit) and vertex to calculate asymmetry score
};

} // namespace lar_content

#endif // #ifndef LAR_GLOBAL_ASYMMETRY_FEATURE_TOOL_H
