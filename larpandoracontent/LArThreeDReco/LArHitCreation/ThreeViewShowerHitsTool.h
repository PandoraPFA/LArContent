/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h
 *
 *  @brief  Header file for the three view shower hits tool.
 *
 *  $Log: $
 */
#ifndef THREE_VIEW_SHOWER_HITS_TOOL_H
#define THREE_VIEW_SHOWER_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ThreeViewShowerHitsTool class
 */
class ThreeViewShowerHitsTool : public ShowerHitsBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeViewShowerHitsTool();

private:
    void GetShowerHit3D(const pandora::CaloHitVector &caloHitVector1, const pandora::CaloHitVector &caloHitVector2, ProtoHit &protoHit) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_zTolerance; ///< The z tolerance to use when looking for associated calo hits between views
};

} // namespace lar_content

#endif // #ifndef THREE_VIEW_SHOWER_HITS_TOOL_H
