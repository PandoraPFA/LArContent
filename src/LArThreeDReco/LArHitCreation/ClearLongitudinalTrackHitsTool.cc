/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.cc
 * 
 *  @brief  Implementation of the longitudinal track hit creation tool.
 * 
 *  $Log: $
 */


#include "LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h"

using namespace pandora;

namespace lar
{
 
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearLongitudinalTrackHitsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  return LongitudinalTrackHitsBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
