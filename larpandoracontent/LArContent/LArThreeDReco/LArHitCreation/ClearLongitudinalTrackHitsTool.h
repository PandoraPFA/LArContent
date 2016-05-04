/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h
 *
 *  @brief  Header file for the clear longitudinal track hit creation tool.
 *
 *  $Log: $
 */
#ifndef CLEAR_LONGITUDINAL_TRACK_HITS_TOOL_H
#define CLEAR_LONGITUDINAL_TRACK_HITS_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ClearLongitudinalTrackHitsTool class
 */
class ClearLongitudinalTrackHitsTool : public LongitudinalTrackHitsBaseTool
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
        const pandora::CartesianVector &vtx3D, const pandora::CartesianVector &end3D, pandora::CartesianVector &position3D, float &chiSquared) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearLongitudinalTrackHitsTool::Factory::CreateAlgorithmTool() const
{
    return new ClearLongitudinalTrackHitsTool();
}

} // namespace lar_content

#endif // #ifndef CLEAR_LONGITUDINAL_TRACK_HITS_TOOL_H
