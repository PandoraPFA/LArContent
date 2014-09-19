/**
 *  @file   LArContent/include/LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.h
 *
 *  @brief  Header file for the clear tracks tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_LONGITUDINAL_TRACKS_TOOL_H
#define CLEAR_LONGITUDINAL_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearLongitudinalTracksTool class
 */
class ClearLongitudinalTracksTool : public LongitudinalTensorTool
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
    ClearLongitudinalTracksTool();

    bool Run(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearLongitudinalTracksTool::Factory::CreateAlgorithmTool() const
{
    return new ClearLongitudinalTracksTool();
}

} // namespace lar_content

#endif // #ifndef CLEAR_LONGITUDINAL_TRACKS_TOOL_H
