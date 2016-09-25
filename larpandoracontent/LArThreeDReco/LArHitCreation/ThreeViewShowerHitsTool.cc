/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.cc
 * 
 *  @brief  Implementation of the three view shower hits tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h"

using namespace pandora;

namespace lar_content
{

ThreeViewShowerHitsTool::ThreeViewShowerHitsTool() :
    m_zTolerance(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowerHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const CaloHitVector &caloHitVector1, const CaloHitVector &caloHitVector2,
    CartesianVector &position3D, float &chiSquared) const
{
    if (caloHitVector1.empty() || caloHitVector2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType hitType1(caloHitVector1.at(0)->GetHitType());
    const HitType hitType2(caloHitVector2.at(0)->GetHitType());

    const HitType hitType2D(pCaloHit2D->GetHitType());
    const float position2D(pCaloHit2D->GetPositionVector().GetZ());

    bool positionFound(false);

    for (const CaloHit *const pCaloHit1 : caloHitVector1)
    {
        const CartesianVector &position1(pCaloHit1->GetPositionVector());
        const float prediction(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2D, hitType1, position2D, position1.GetZ()));

        for (const CaloHit *const pCaloHit2 : caloHitVector2)
        {
            const CartesianVector &position2(pCaloHit2->GetPositionVector());

            if (std::fabs(position2.GetZ() - prediction) > m_zTolerance)
                continue;

            CartesianVector thisPosition3D(0.f, 0.f, 0.f);
            float thisChiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, hitType1, hitType2, position1, position2, thisPosition3D, thisChiSquared); 

            if (thisChiSquared < chiSquared)
            {
                chiSquared = thisChiSquared;
                position3D = thisPosition3D;
                positionFound = true;
            }
        }
    }

    if (!positionFound)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewShowerHitsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ZTolerance", m_zTolerance));

    return ShowerHitsBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
