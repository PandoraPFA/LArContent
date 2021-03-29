/**
 *  @file   larpandoracontent/LArThreeDReco/LArTrackMatching/TwoViewSimpleTracksTool.h
 *
 *  @brief  Header file for the simple showers tool class.
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_SIMPLE_TRACKS_TOOL_H
#define TWO_VIEW_SIMPLE_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewSimpleTracksTool class
 */
class TwoViewSimpleTracksTool : public TransverseMatrixTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewSimpleTracksTool();

    bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

private:
    /**
     *  @brief  Find best track match as a simple way to (try to) resolve ambiguities in the matrix
     *
     *  @param  overlapMatrix the overlap matrix
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindBestTrack(const MatrixType &overlapMatrix, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Whether a provided (iterator to a) matrix element passes the selection cuts for particle creation
     *
     *  @param  eIter the iterator to the matrix element
     */
    bool PassesElementCuts(MatrixType::ElementList::const_reverse_iterator eIter) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minMatchedFraction;              ///< The min matched sampling point fraction for particle creation
    float m_minMatchingScore;                ///< The min global matching score for particle creation
    unsigned int m_minMatchedSamplingPoints; ///< The min number of matched sampling points for particle creation
    float m_minXOverlapFraction;             ///< The min x overlap fraction (in each view) for particle creation
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_SIMPLE_TRACKS_TOOL_H
