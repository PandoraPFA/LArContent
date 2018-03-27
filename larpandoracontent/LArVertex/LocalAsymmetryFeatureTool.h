/**
 *  @file   larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h
 *
 *  @brief  Header file for the local asymmetry feature tool class.
 *
 *  $Log: $
 */
#ifndef LAR_LOCAL_ASYMMETRY_FEATURE_TOOL_H
#define LAR_LOCAL_ASYMMETRY_FEATURE_TOOL_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  LocalAsymmetryFeatureTool class
 */
class LocalAsymmetryFeatureTool : public VertexSelectionBaseAlgorithm::VertexFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    LocalAsymmetryFeatureTool();

    /**
     *  @brief  Run the tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pVertex address of the vertex
     *  @param  slidingFitDataListMap map of the sliding fit data lists
     *
     *  @return the energy kick feature
     */
    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm * const pAlgorithm, const pandora::Vertex * const pVertex,
        const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap,const VertexSelectionBaseAlgorithm::ClusterListMap &,
        const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the local asymmetry feature in a given view
     *
     *  @param  vertexPosition2D the vertex position projected into this view
     *  @param  slidingFitDataList the list of sliding fit data objects in this view
     *
     *  @return the local asymmetry feature
     */
    float GetLocalAsymmetryForView(const pandora::CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList) const;

    /**
     *  @brief  Increment the asymmetry parameters
     *
     *  @param  weight the weight for the increment
     *  @param  clusterDirection the cluster direction
     *  @param  localWeightedDirectionSum the current value of the vector
     *
     *  @return whether the energy asymmetry score is still viable
     */
    bool IncrementAsymmetryParameters(const float weight, const pandora::CartesianVector &clusterDirection,
        pandora::CartesianVector &localWeightedDirectionSum) const;

    /**
     *  @brief  Calculate the local asymmetry feature
     *
     *  @param  useEnergyMetrics whether to use energy-based rather than hit-counting-based metrics
     *  @param  vertexPosition2D the vertex position
     *  @param  asymmetryClusters the clusters to use to calculate the asymmetry
     *  @param  localWeightedDirectionSum the local event axis
     *
     *  @return the local asymmetry feature
     */
    float CalculateLocalAsymmetry(const bool useEnergyMetrics, const pandora::CartesianVector &vertexPosition2D,
        const pandora::ClusterVector &asymmetryClusters, const pandora::CartesianVector &localWeightedDirectionSum) const;

    float           m_maxAsymmetryDistance;     ///< The max distance between cluster (any hit) and vertex to calculate asymmetry score
    float           m_minAsymmetryCosAngle;     ///< The min opening angle cosine used to determine viability of asymmetry score
    unsigned int    m_maxAsymmetryNClusters;    ///< The max number of associated clusters to calculate the asymmetry
};

} // namespace lar_content

#endif // #ifndef LAR_LOCAL_ASYMMETRY_FEATURE_TOOL_H
