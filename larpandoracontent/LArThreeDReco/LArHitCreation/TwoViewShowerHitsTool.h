/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h
 *
 *  @brief  Header file for the two view shower hits tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_SHOWER_HITS_TOOL_H
#define TWO_VIEW_SHOWER_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  TwoViewShowerHitsTool class
 */
class TwoViewShowerHitsTool : public ShowerHitsBaseTool
{
private:
    void GetShowerHit3D(const pandora::CaloHitVector &caloHitVector1, const pandora::CaloHitVector &caloHitVector2, ProtoHit &protoHit) const;

    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a list of candidate matched
     *          hits in one of the other two views
     *
     *  @param  caloHitVector the vector of candidate hits in another view
     *  @param  protoHit to receive the populated proto hit
     */
    void GetShowerHit3D(const pandora::CaloHitVector &caloHitVector, ProtoHit &protoHit) const;
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_SHOWER_HITS_TOOL_H
