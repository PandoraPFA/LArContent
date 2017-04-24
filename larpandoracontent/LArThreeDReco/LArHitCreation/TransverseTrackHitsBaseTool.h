/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h
 * 
 *  @brief  Header file for the transverse track hits base tool.
 * 
 *  $Log: $
 */
#ifndef TRANSVERSE_TRACK_HITS_BASE_TOOL_H
#define TRANSVERSE_TRACK_HITS_BASE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

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
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    virtual void CreateThreeDHits(const pandora::CaloHitVector &inputTwoDHits, const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and sliding linear fits in the other two views
     * 
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetThreeDPosition(const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHit &protoHit) const = 0;

    /**
     *  @brief  Calculate an additional contribution to the chi-squared based on the steepness of the track
     * 
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  protoHit to receive the modified proto hit
     */
    virtual void GetTransverseChi2(const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHit &protoHit) const;
        
    /**
     *  @brief  Calculate an additional contribution to the chi-squared based on the steepness of the track
     * 
     *  @param  position2D the calculated two dimensional position
     *  @param  fitResult the sliding fit to the track
     */
    virtual double GetTransverseChi2(const pandora::CartesianVector &position2D, const TwoDSlidingFitResult &fitResult) const;
};

} // namespace lar_content

#endif // #ifndef TRANSVERSE_TRACK_HITS_BASE_TOOL_H
