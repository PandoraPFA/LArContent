/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h
 *
 *  @brief  Header file for the clear longitudinal track hit creation tool.
 *
 *  $Log: $
 */
#ifndef CLEAR_LONGITUDINAL_TRACK_HITS_TOOL_H
#define CLEAR_LONGITUDINAL_TRACK_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ClearLongitudinalTrackHitsTool class
 */
class ClearLongitudinalTrackHitsTool : public LongitudinalTrackHitsBaseTool
{
private:
    void GetLongitudinalTrackHit3D(const MatchedSlidingFitMap &matchedSlidingFitMap, const pandora::CartesianVector &vtx3D,
        const pandora::CartesianVector &end3D, ProtoHit &protoHit) const;
};

} // namespace lar_content

#endif // #ifndef CLEAR_LONGITUDINAL_TRACK_HITS_TOOL_H
