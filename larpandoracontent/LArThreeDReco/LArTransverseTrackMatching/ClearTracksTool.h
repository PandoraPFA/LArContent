/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h
 *
 *  @brief  Header file for the clear tracks tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_TRACKS_TOOL_H
#define CLEAR_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearTracksTool class
 */
class ClearTracksTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ClearTracksTool();

    bool Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;

    float m_minMatchedFraction;  ///< The min matched sampling point fraction for particle creation
    float m_minXOverlapFraction; ///< The min x overlap fraction (in each view) for particle creation
    bool m_visualize;            ///< Visualize the clear track matches
};

} // namespace lar_content

#endif // #ifndef CLEAR_TRACKS_TOOL_H
