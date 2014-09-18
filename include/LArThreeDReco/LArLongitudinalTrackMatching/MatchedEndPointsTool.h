/**
 *  @file   LArContent/include/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h
 *
 *  @brief  Header file for the matched end points tool class.
 *
 *  $Log: $
 */
#ifndef MATCHED_END_POINTS_TOOL_H
#define MATCHED_END_POINTS_TOOL_H 1

#include "LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  MatchedEndPointsTool class
 */
class MatchedEndPointsTool : public LongitudinalTensorTool
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
    MatchedEndPointsTool();

    bool Run(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find matched tracks, hidden by ambiguities in the tensor
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindMatchedTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    float           m_maxEndPointChi2;                    ///< The max chi2 of matched vertex and end points for particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *MatchedEndPointsTool::Factory::CreateAlgorithmTool() const
{
    return new MatchedEndPointsTool();
}

} // namespace lar_content

#endif // #ifndef MATCHED_END_POINTS_TOOL_H
