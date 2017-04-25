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
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

private:
    void GetLongitudinalTrackHit3D(const MatchedSlidingFitMap &matchedSlidingFitMap, const pandora::CartesianVector &vtx3D,
        const pandora::CartesianVector &end3D, ProtoHit &protoHit) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *MultiValuedLongitudinalTrackHitsTool::Factory::CreateAlgorithmTool() const
{
    return new MultiValuedLongitudinalTrackHitsTool();
}

} // namespace lar_content

#endif // #ifndef MULTI_VALUED_LONGITUDINAL_TRACK_HITS_TOOL_H
