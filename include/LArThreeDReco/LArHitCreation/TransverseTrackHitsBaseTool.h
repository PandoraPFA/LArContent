/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h
 * 
 *  @brief  Header file for the transverse track hits base tool.
 * 
 *  $Log: $
 */
#ifndef TRANSVERSE_TRACK_HITS_BASE_TOOL_H
#define TRANSVERSE_TRACK_HITS_BASE_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  TransverseTrackHitsBaseTool class
 */
class TransverseTrackHitsBaseTool : public TrackHitsBaseTool
{
protected:
    /**
     *  @brief  Create three dimensional hits, using an input list of two dimensional hits and two associated sliding fit results
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    virtual void CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::CaloHitList &inputTwoDHits, 
        const MatchedSlidingFitMap &matchedSlidingFitMap, pandora::CaloHitList &newThreeDHits) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and sliding linear fits in the other two views
     * 
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetThreeDPosition(const pandora::CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
        pandora::CartesianVector &position3D, float &chiSquared) const = 0;

    /**
     *  @brief  Calculate an additional contribution to the chi-squared based on the steepness of the track
     * 
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  position3D the calculated three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetTransverseChi2(const pandora::CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
        const pandora::CartesianVector &position3D, float &chiSquared) const;
        
    /**
     *  @brief  Calculate an additional contribution to the chi-squared based on the steepness of the track
     * 
     *  @param  position2D the calculated two dimensional position
     *  @param  fitResult the sliding fit to the track
     */
    virtual float GetTransverseChi2(const pandora::CartesianVector &position2D, const TwoDSlidingFitResult &fitResult) const;
};

} // namespace lar_content

#endif // #ifndef TRANSVERSE_TRACK_HITS_BASE_TOOL_H
