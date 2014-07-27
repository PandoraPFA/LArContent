/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.cc
 * 
 *  @brief  Implementation of the longitudinal track hit creation tool.
 * 
 *  $Log: $
 */

#include "LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h"

using namespace pandora;

namespace lar
{
 
void LongitudinalTrackHitsBaseTool::CreateThreeDHits(ThreeDHitCreationAlgorithm */*pAlgorithm*/, const CaloHitList &/*inputTwoDHits*/, 
    const MatchedSlidingFitMap &/*matchedSlidingFitMap*/, CaloHitList &/*newThreeDHits*/, CaloHitList &/*omittedTwoDHits*/) const
{   
  
  // TODO --- 
 
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalTrackHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return TrackHitsBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
