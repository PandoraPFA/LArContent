/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h
 * 
 *  @brief  Header file for the longitudinal track hit creation tool.
 * 
 *  $Log: $
 */
#ifndef LONGITUDINAL_TRACK_HITS_BASE_TOOL_H
#define LONGITUDINAL_TRACK_HITS_BASE_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

namespace lar
{

/**
 *  @brief  LongitudinalTrackHitsBaseTool class
 */
class LongitudinalTrackHitsBaseTool : public TrackHitsBaseTool
{
protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional hits, using an input list of two dimensional hits and two associated sliding fit results
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  matchedSlidingFitMap the sliding fit results for each view
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    void CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::CaloHitList &inputTwoDHits, 
        const MatchedSlidingFitMap &matchedSlidingFitMap, pandora::CaloHitList &newThreeDHits) const;

 

};

} // namespace lar

#endif // #ifndef LONGITUDINAL_TRACK_HITS_BASE_TOOL_H
