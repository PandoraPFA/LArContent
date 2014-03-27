/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/OvershootTracksTool.h
 * 
 *  @brief  Header file for the overshoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef OVERSHOOT_TRACKS_TOOL_H
#define OVERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  OvershootTracksTool class
 */
class OvershootTracksTool : public TensorManipulationTool
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
    pandora::StatusCode Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find overshoot tracks
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindOvershootTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *OvershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new OvershootTracksTool();
}

} // namespace lar

#endif // #ifndef OVERSHOOT_TRACKS_TOOL_H
