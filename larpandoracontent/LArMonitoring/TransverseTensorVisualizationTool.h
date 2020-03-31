/**
 *  @file   larpandoracontent/LArMonitoring/TransverseTensorVisualizationTool.h
 *
 *  @brief  Header file for the transverse tensor visualization tool class.
 *
 *  $Log: $
 */
#ifndef TRANSVERSE_TENSOR_VISUALIZATION_TOOL_H
#define TRANSVERSE_TENSOR_VISUALIZATION_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TransverseTensorVisualizationTool class
 */
class TransverseTensorVisualizationTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TransverseTensorVisualizationTool();

    bool Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minClusterConnections;        ///< The minimum number of cluster connections for display
    bool            m_ignoreUnavailableClusters;    ///< Whether to ignore (skip-over) unavailable clusters in the tensor
    bool            m_showEachIndividualElement;    ///< Whether to draw each individual tensor element
    bool            m_showContext;                  ///< Whether to show input cluster lists to add context to tensor elements
};

} // namespace lar_content

#endif // #ifndef TRANSVERSE_TENSOR_VISUALIZATION_TOOL_H
