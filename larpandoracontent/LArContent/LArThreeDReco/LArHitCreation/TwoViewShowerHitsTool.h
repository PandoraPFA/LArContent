/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h
 *
 *  @brief  Header file for the two view shower hits tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_SHOWER_HITS_TOOL_H
#define TWO_VIEW_SHOWER_HITS_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  TwoViewShowerHitsTool class
 */
class TwoViewShowerHitsTool : public ShowerHitsBaseTool
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

    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a list of candidate matched
     *          hits in one of the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  caloHitList the list of candidate hits in another view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetThreeDPosition(const pandora::CaloHit *const pCaloHit2D, const pandora::CaloHitList &caloHitList, pandora::CartesianVector &position3D,
        float &chiSquared) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TwoViewShowerHitsTool::Factory::CreateAlgorithmTool() const
{
    return new TwoViewShowerHitsTool();
}

} // namespace lar_content

#endif // #ifndef TWO_VIEW_SHOWER_HITS_TOOL_H
