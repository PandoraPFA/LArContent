/**
 *  @file   larpandoracontent/LArMonitoring/LongitudinalTensorVisualizationTool.h
 *
 *  @brief  Header file for the longitudinal tensor visualization tool class.
 *
 *  $Log: $
 */
#ifndef LONGITUDINAL_TENSOR_VISUALIZATION_TOOL_H
#define LONGITUDINAL_TENSOR_VISUALIZATION_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeViewLongitudinalTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  LongitudinalTensorVisualizationTool class
 */
class LongitudinalTensorVisualizationTool : public LongitudinalTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    LongitudinalTensorVisualizationTool();

    bool Run(ThreeViewLongitudinalTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_minClusterConnections; ///< The minimum number of cluster connections for display
    bool m_ignoreUnavailableClusters;     ///< Whether to ignore (skip-over) unavailable clusters in the tensor
    bool m_showEachIndividualElement;     ///< Whether to draw each individual tensor element
    bool m_showContext;                   ///< Whether to show input cluster lists to add context to tensor elements
};

} // namespace lar_content

#endif // #ifndef LONGITUDINAL_TENSOR_VISUALIZATION_TOOL_H
