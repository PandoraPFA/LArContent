/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h
 * 
 *  @brief  Header file for the multivalued transverse track hit creation tool.
 * 
 *  $Log: $
 */
#ifndef MULTI_VALUED_TRANSVERSE_TRACK_HITS_TOOL_H
#define MULTI_VALUED_TRANSVERSE_TRACK_HITS_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h"

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
    void GetThreeDPosition(const pandora::CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
        pandora::CartesianVector &position3D, float &chiSquared) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *MultiValuedTransverseTrackHitsTool::Factory::CreateAlgorithmTool() const
{
    return new MultiValuedTransverseTrackHitsTool();
}

} // namespace lar_content

#endif // #ifndef MULTI_VALUED_TRANSVERSE_TRACK_HITS_TOOL_H
