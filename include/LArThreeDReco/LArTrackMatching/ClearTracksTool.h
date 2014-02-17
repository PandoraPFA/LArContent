/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ClearTracksTool.h
 * 
 *  @brief  Header file for the clear tracks tool class.
 * 
 *  $Log: $
 */
#ifndef CLEAR_TRACKS_TOOL_H
#define CLEAR_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  ClearTracksTool class
 */
class ClearTracksTool : public TensorManipulationTool
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
    pandora::StatusCode Run(OverlapTensor<TrackOverlapResult> &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearTracksTool::Factory::CreateAlgorithmTool() const
{
    return new ClearTracksTool();
}

} // namespace lar

#endif // #ifndef CLEAR_TRACKS_TOOL_H
