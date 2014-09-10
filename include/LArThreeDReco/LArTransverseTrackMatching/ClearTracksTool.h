/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h
 * 
 *  @brief  Header file for the clear tracks tool class.
 * 
 *  $Log: $
 */
#ifndef CLEAR_TRACKS_TOOL_H
#define CLEAR_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearTracksTool class
 */
class ClearTracksTool : public TransverseTensorTool
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

    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     * 
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearTracksTool::Factory::CreateAlgorithmTool() const
{
    return new ClearTracksTool();
}

} // namespace lar_content

#endif // #ifndef CLEAR_TRACKS_TOOL_H
