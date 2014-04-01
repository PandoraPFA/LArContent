/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/MissingTracksTool.h
 * 
 *  @brief  Header file for the missing tracks tool class.
 * 
 *  $Log: $
 */
#ifndef MISSING_TRACKS_TOOL_H
#define MISSING_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  MissingTracksTool class
 */
class MissingTracksTool : public TensorManipulationTool
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

private:
    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *MissingTracksTool::Factory::CreateAlgorithmTool() const
{
    return new MissingTracksTool();
}

} // namespace lar

#endif // #ifndef MISSING_TRACKS_TOOL_H
