/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewClearTracksTool.h
 *
 *  @brief  Header file for the two view clear tracks tool class.
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_CLEAR_TRACKS_TOOL_H
#define TWO_VIEW_CLEAR_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewClearTracksTool class
 */
class TwoViewClearTracksTool : public TransverseMatrixTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewClearTracksTool();

    bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList, bool &particlesMade) const;

    float m_minXOverlapFraction;       ///< The min x overlap fraction value for particle creation
    float m_minMatchingScore;          ///< The min global matching score for particle creation
    float m_minLocallyMatchedFraction; ///< The min locally matched fraction for particle creation
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_CLEAR_TRACKS_TOOL_H
