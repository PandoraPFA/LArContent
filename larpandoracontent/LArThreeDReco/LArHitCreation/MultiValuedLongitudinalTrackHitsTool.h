/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.h
 *
 *  @brief  Header file for the multivalued longitudinal track hit creation tool.
 *
 *  $Log: $
 */
#ifndef MULTI_VALUED_LONGITUDINAL_TRACK_HITS_TOOL_H
#define MULTI_VALUED_LONGITUDINAL_TRACK_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  MultiValuedLongitudinalTrackHitsTool class
 */
class MultiValuedLongitudinalTrackHitsTool : public LongitudinalTrackHitsBaseTool
{
private:
    void GetLongitudinalTrackHit3D(const MatchedSlidingFitMap &matchedSlidingFitMap, const pandora::CartesianVector &vtx3D,
        const pandora::CartesianVector &end3D, ProtoHit &protoHit) const;
};

} // namespace lar_content

#endif // #ifndef MULTI_VALUED_LONGITUDINAL_TRACK_HITS_TOOL_H
