/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/TensorManipulationTool.h
 * 
 *  @brief  Header file for the tensor manipulation tool class.
 * 
 *  $Log: $
 */
#ifndef TENSOR_MANIPULATION_TOOL_H
#define TENSOR_MANIPULATION_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "LArObjects/LArOverlapTensor.h"

namespace lar
{

/**
 *  @brief  TensorManipulationTool class
 */
class TensorManipulationTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     * 
     *  @param  overlapTensor the overlap tensor
     */
    virtual pandora::StatusCode Run(OverlapTensor<TrackOverlapResult> &overlapTensor) = 0;
};

} // namespace lar

#endif // #ifndef TENSOR_MANIPULATION_TOOL_H
