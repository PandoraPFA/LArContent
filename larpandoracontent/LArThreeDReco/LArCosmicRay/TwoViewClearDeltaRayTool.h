/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewClearDeltaRayTool.h
 *
 *  @brief  Header file for the two view clear delta ray tool class
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_CLEAR_DELTA_RAY_TOOL_H
#define TWO_VIEW_CLEAR_DELTA_RAY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewClearDeltaRayTool class
 */
class TwoViewClearDeltaRayTool : public DeltaRayMatrixTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewClearDeltaRayTool();

    bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList, bool &particlesMade) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_CLEAR_TRACKS_TOOL_H
