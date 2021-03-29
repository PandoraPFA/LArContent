/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h
 *
 *  @brief  Header file for the matched end points tool class.
 *
 *  $Log: $
 */
#ifndef MATCHED_END_POINTS_TOOL_H
#define MATCHED_END_POINTS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeViewLongitudinalTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  MatchedEndPointsTool class
 */
class MatchedEndPointsTool : public LongitudinalTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    MatchedEndPointsTool();

    bool Run(ThreeViewLongitudinalTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find matched tracks, hidden by ambiguities in the tensor
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindMatchedTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Sort tensor elements by chi-squared
     *
     *  @param  lhs the first tensor element
     *  @param  rhs the second tensor element
     *
     *  @return boolean
     */
    static bool SortByChiSquared(const TensorType::Element &lhs, const TensorType::Element &rhs);

    float m_minMatchedFraction; ///< The min matched sampling point fraction for particle creation
    float m_maxEndPointChi2;    ///< The max chi2 of matched vertex and end points for particle creation
};

} // namespace lar_content

#endif // #ifndef MATCHED_END_POINTS_TOOL_H
