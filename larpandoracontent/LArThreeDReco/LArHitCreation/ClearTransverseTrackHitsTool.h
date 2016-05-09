/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h
 *
 *  @brief  Header file for the clear transverse track hit creation tool.
 *
 *  $Log: $
 */
#ifndef CLEAR_TRANSVERSE_TRACK_HITS_TOOL_H
#define CLEAR_TRANSVERSE_TRACK_HITS_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ClearTransverseTrackHitsTool class
 */
class ClearTransverseTrackHitsTool : public TransverseTrackHitsBaseTool
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

inline pandora::AlgorithmTool *ClearTransverseTrackHitsTool::Factory::CreateAlgorithmTool() const
{
    return new ClearTransverseTrackHitsTool();
}

} // namespace lar_content

#endif // #ifndef CLEAR_TRANSVERSE_TRACK_HITS_TOOL_H
