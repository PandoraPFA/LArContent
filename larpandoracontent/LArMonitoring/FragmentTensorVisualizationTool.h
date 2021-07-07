/**
 *  @file   larpandoracontent/LArMonitoring/FragmentTensorVisualizationTool.h
 *
 *  @brief  Header file for the fragment tensor visualization tool class.
 *
 *  $Log: $
 */
#ifndef FRAGMENT_TENSOR_VISUALIZATION_TOOL_H
#define FRAGMENT_TENSOR_VISUALIZATION_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTrackFragments/ThreeViewTrackFragmentsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  FragmentTensorVisualizationTool class
 */
class FragmentTensorVisualizationTool : public FragmentTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    FragmentTensorVisualizationTool();

    bool Run(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_minClusterConnections; ///< The minimum number of cluster connections for display
    bool m_ignoreUnavailableClusters;     ///< Whether to ignore (skip-over) unavailable clusters in the tensor
    bool m_showEachIndividualElement;     ///< Whether to draw each individual tensor element
    bool m_showContext;                   ///< Whether to show input cluster lists to add context to tensor elements
};

} // namespace lar_content

#endif // #ifndef FRAGMENT_TENSOR_VISUALIZATION_TOOL_H
