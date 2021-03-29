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
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and sliding linear fits in the other two views
     *
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetTransverseTrackHit3D(const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHit &protoHit) const = 0;

    virtual void GetTrackHits3D(
        const pandora::CaloHitVector &inputTwoDHits, const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Calculate an additional contribution to the chi-squared based on the steepness of the track
     *
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  protoHit to receive the modified proto hit
     */
    virtual void AddTransverseChi2(const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHit &protoHit) const;

private:
    /**
     *  @brief  Calculate an additional contribution to the chi-squared based on the steepness of the track
     *
     *  @param  position2D the calculated two dimensional position
     *  @param  fitResult the sliding fit to the track
     */
    double GetTransverseChi2(const pandora::CartesianVector &position2D, const TwoDSlidingFitResult &fitResult) const;
};

} // namespace lar_content

#endif // #ifndef TRANSVERSE_TRACK_HITS_BASE_TOOL_H
