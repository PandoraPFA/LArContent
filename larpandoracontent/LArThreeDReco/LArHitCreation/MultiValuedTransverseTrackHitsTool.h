/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h
 * 
 *  @brief  Header file for the multivalued transverse track hit creation tool.
 * 
 *  $Log: $
 */
#ifndef MULTI_VALUED_TRANSVERSE_TRACK_HITS_TOOL_H
#define MULTI_VALUED_TRANSVERSE_TRACK_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  MultiValuedTransverseTrackHitsTool class
 */
class MultiValuedTransverseTrackHitsTool : public TransverseTrackHitsBaseTool
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
    void GetTransverseTrackHit3D(const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHit &protoHit) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *MultiValuedTransverseTrackHitsTool::Factory::CreateAlgorithmTool() const
{
    return new MultiValuedTransverseTrackHitsTool();
}

} // namespace lar_content

#endif // #ifndef MULTI_VALUED_TRANSVERSE_TRACK_HITS_TOOL_H
