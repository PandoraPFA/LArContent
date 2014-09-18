/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.h
 * 
 *  @brief  Header file for the shower tensor visualization tool class.
 * 
 *  $Log: $
 */
#ifndef SHOWER_TENSOR_VISUALIZATION_TOOL_H
#define SHOWER_TENSOR_VISUALIZATION_TOOL_H 1

#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerTensorVisualizationTool class
 */
class ShowerTensorVisualizationTool : public ShowerTensorTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    ShowerTensorVisualizationTool();

    bool Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minClusterConnections;        ///< The minimum number of cluster connections for display
    bool            m_ignoreUnavailableClusters;    ///< Whether to ignore (skip-over) unavailable clusters in the tensor
    bool            m_showEachIndividualElement;    ///< Whether to draw each individual tensor element
    bool            m_showContext;                  ///< Whether to show input cluster lists to add context to tensor elements
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ShowerTensorVisualizationTool::Factory::CreateAlgorithmTool() const
{
    return new ShowerTensorVisualizationTool();
}

} // namespace lar_content

#endif // #ifndef SHOWER_TENSOR_VISUALIZATION_TOOL_H
