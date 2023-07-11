/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h
 *
 *  @brief  Header file for the peak direction finder tool class.
  *
 *  $Log: $
 */
#ifndef LAR_SHOWER_SPINE_FINDER_TOOL_H
#define LAR_SHOWER_SPINE_FINDER_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

class ShowerSpineFinderTool : public pandora::AlgorithmTool
{
public:
    ShowerSpineFinderTool();

    pandora::StatusCode Run(const pandora::CartesianVector &nuVertex3D, const pandora::CaloHitList *const pViewHitList, const pandora::HitType hitType,
        const pandora::CartesianVector &peakDirection, pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Perform a running fit to collect the hits of the shower spine
     *
     *  @param  pViewHitList the 2D event hit list
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  initialDirection the initial direction of the pathway
     *  @param  unavailableHitList protected hits that cannot be collected
     *  @param  showerSpineHitList the output list of shower spine hits
     */
    void FindShowerSpine(const pandora::CaloHitList *const pViewHitList, const pandora::CartesianVector &nuVertex2D,
        const pandora::CartesianVector &initialDirection, pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList) const;

    /**
     *  @brief  Perform a running fit step: collect hits which lie close to the shower spine projection
     *
     *  @param  extrapolatedFit the fit to the hitherto collected hits
     *  @param  extrapolatedStartPosition the shower spine projection start position
     *  @param  extrapolatedEndPosition the shower spine projection end position
     *  @param  extrapolatedDirection the shower spine projection direction
     *  @param  isEndDownstream whether the shower direction is downstream (in Z) of the neutrino vertex
     *  @param  pViewHitList the 2D event hit list
     *  @param  runningFitPositionVector the vector of the hitherto collected hit positions
     *  @param  unavailableHitList protected hits that cannot be collected
     *  @param  showerSpineHitList the output list of shower spine hits
     *
     *  @return whether any hits were collected in the running fit step
     */
    bool CollectSubsectionHits(const TwoDSlidingFitResult &extrapolatedFit, const pandora::CartesianVector &extrapolatedStartPosition,
        const pandora::CartesianVector &extrapolatedEndPosition, const pandora::CartesianVector &extrapolatedDirection,
        const bool isEndDownstream, const pandora::CaloHitList *const pViewHitList, pandora::CartesianPointVector &runningFitPositionVector,
        pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList) const;

    /**
     *  @brief  Determine whether a hit lies close to the shower spine projection
     *
     *  @param  hitPosition the hit position
     *  @param  lineStart the shower spine projection start position
     *  @param  lineDirection the shower spine projection direction
     *  @param  distanceToLine the comparison distance for 'is close'
     *
     *  @return whether the hit is close to the shower spine projection
     */
    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart,
        const pandora::CartesianVector &lineDirection, const float distanceToLine) const;

    /**
     *  @brief  Add to the shower spine the connecting hits
     *
     *  @param  collectedHits the input list of close hits
     *  @param  extrapolatedStartPosition the shower spine projection start position
     *  @param  extrapolatedDirection the shower spine projection direction
     *  @param  runningFitPositionVector the vector of the collected hit positions
     *  @param  showerSpineHitList the list of collected shower spine hits
     */
    void CollectConnectedHits(const pandora::CaloHitList &collectedHits, const pandora::CartesianVector &extrapolatedStartPosition,
        const pandora::CartesianVector &extrapolatedDirection, pandora::CartesianPointVector &runningFitPositionVector,
        pandora::CaloHitList &showerSpineHitList) const;

    /**
     *  @brief  Find the smallest distance between a position and a list of other positions
     *
     *  @param  position the input position
     *  @param  testPositions the list of other positions
     *
     *  @return the closest distance
     */
    float GetClosestDistance(const pandora::CartesianVector &position, const pandora::CartesianPointVector &testPositions) const;

    unsigned int m_hitThresholdForSpine;           ///< The hit threshold for a significant spine
    float m_growingFitInitialLength;               ///< The first step distance
    float m_initialFitDistanceToLine;              ///< The max. proximity to the spine projection for collection in the first step
    unsigned int m_minInitialHitsFound;            ///< The min. number of hits collected in the first step for continuation
    unsigned int m_maxFittingHits;                 ///< The number of hits to consider in the running fit
    unsigned int m_localSlidingFitWindow;          ///< The standard sliding fit window for spine fits
    float m_growingFitSegmentLength;               ///< The standard step distance
    unsigned int m_highResolutionSlidingFitWindow; ///< The high resolution sliding fit window for spine fits
    float m_distanceToLine;                        ///< The max. proximity to the spine projection for collection
    float m_hitConnectionDistance;                 ///< The max. separation between connected hits
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_SHOWER_SPINE_FINDER_TOOL_H
