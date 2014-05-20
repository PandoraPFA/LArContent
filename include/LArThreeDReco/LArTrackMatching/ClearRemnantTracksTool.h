/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ClearRemnantTracksTool.h
 *
 *  @brief  Header file for the clear tracks tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_REMNANT_TRACKS_TOOL_H
#define CLEAR_REMNANT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  ClearRemnantTracksTool class
 */
class ClearRemnantTracksTool : public RemnantTensorTool
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

    bool Run(ThreeDRemnantTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeDRemnantTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearRemnantTracksTool::Factory::CreateAlgorithmTool() const
{
    return new ClearRemnantTracksTool();
}

} // namespace lar

#endif // #ifndef CLEAR_REMNANT_TRACKS_TOOL_H
