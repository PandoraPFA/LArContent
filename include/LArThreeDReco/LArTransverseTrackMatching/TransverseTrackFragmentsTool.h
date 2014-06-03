/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/TransverseTrackFragmentsTool.h
 *
 *  @brief  Header file for the transverse track fragments tool class.
 *
 *  $Log: $
 */
#ifndef TRANSVERSE_TRACK_FRAGMENTS_TOOL_H
#define TRANSVERSE_TRACK_FRAGMENTS_TOOL_H 1

#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTrackFragmentsAlg.h"

namespace lar
{

/**
 *  @brief  TransverseTrackFragmentsTool class
 */
class TransverseTrackFragmentsTool : public TransverseFragmentTensorTool
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

    bool Run(ThreeDTransverseTrackFragmentsAlg *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TransverseTrackFragmentsTool::Factory::CreateAlgorithmTool() const
{
    return new TransverseTrackFragmentsTool();
}

} // namespace lar

#endif // #ifndef TRANSVERSE_TRACK_FRAGMENTS_TOOL_H
