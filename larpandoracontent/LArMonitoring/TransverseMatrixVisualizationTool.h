/**
 *  @file   larpandoracontent/LArMonitoring/TransverseMatrixVisualizationTool.h
 *
 *  @brief  Header file for the transverse matrix visualization tool class.
 *
 *  $Log: $
 */
#ifndef TRANSVERSE_MATRIX_VISUALIZATION_TOOL_H
#define TRANSVERSE_MATRIX_VISUALIZATION_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TransverseMatrixVisualizationTool class
 */
class TransverseMatrixVisualizationTool : public TransverseMatrixTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TransverseMatrixVisualizationTool();

    bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minClusterConnections;                   ///< The minimum number of cluster connections for display
    bool            m_ignoreUnavailableClusters;               ///< Whether to ignore (skip-over) unavailable clusters in the matrix
    bool            m_showEachIndividualElement;               ///< Whether to draw each individual matrix element
    bool            m_showOnlyTrueMatchIndividualElements;     ///< Whether to draw only truly matching individual matrix elements

};

} // namespace lar_content

#endif // #ifndef TRANSVERSE_MATRIX_VISUALIZATION_TOOL_H
