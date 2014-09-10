/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h
 *
 *  @brief  Header file for the three view shower hits tool.
 *
 *  $Log: $
 */
#ifndef THREE_VIEW_SHOWER_HITS_TOOL_H
#define THREE_VIEW_SHOWER_HITS_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ThreeViewShowerHitsTool class
 */
class ThreeViewShowerHitsTool : public ShowerHitsBaseTool
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
    void GetThreeDPosition(const pandora::CaloHit *const pCaloHit2D, const pandora::CaloHitList &caloHitList1, const pandora::CaloHitList &caloHitList2,
        pandora::CartesianVector &position3D, float &chiSquared) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ThreeViewShowerHitsTool::Factory::CreateAlgorithmTool() const
{
    return new ThreeViewShowerHitsTool();
}

} // namespace lar_content

#endif // #ifndef THREE_VIEW_SHOWER_HITS_TOOL_H
