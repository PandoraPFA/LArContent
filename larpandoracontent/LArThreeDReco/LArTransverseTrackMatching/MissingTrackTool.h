/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h
 *
 *  @brief  Header file for the missing track tool class.
 *
 *  $Log: $
 */
#ifndef MISSING_TRACK_TOOL_H
#define MISSING_TRACK_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  MissingTrackTool class
 */
class MissingTrackTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    MissingTrackTool();

    bool Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find missing tracks, due to merging of multiple particle deposits into single hits during hit creation
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindMissingTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    unsigned int    m_minMatchedSamplingPoints;     ///< The min number of matched sampling points for the unavailable tensor element
    float           m_minMatchedFraction;           ///< The min matched sampling point fraction for the unavailable tensor element
    float           m_maxReducedChiSquared;         ///< The max reduced chi squared value for the unavailable tensor element
    float           m_minXOverlapFraction;          ///< The min x overlap fraction for the two available clusters in the tensor element
};

} // namespace lar_content

#endif // #ifndef MISSING_TRACK_TOOL_H
