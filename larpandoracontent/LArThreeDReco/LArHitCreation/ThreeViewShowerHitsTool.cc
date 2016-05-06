/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.cc
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

void ThreeViewShowerHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const CaloHitList &caloHitList1, const CaloHitList &caloHitList2,
    CartesianVector &position3D, float &chiSquared) const
{
    if (caloHitList1.empty() || caloHitList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType hitType1((*caloHitList1.begin())->GetHitType());
    const HitType hitType2((*caloHitList2.begin())->GetHitType());

    const HitType hitType2D(pCaloHit2D->GetHitType());
    const float position2D(pCaloHit2D->GetPositionVector().GetZ());

    bool positionFound(false);

    for (CaloHitList::const_iterator iter1 = caloHitList1.begin(), iterEnd1 = caloHitList1.end(); iter1 != iterEnd1; ++iter1)
    {
        const CartesianVector &position1((*iter1)->GetPositionVector());
        const float prediction(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2D, hitType1, position2D, position1.GetZ()));

        for (CaloHitList::const_iterator iter2 = caloHitList2.begin(), iterEnd2 = caloHitList2.end(); iter2 != iterEnd2; ++iter2)
        {
            const CartesianVector &position2((*iter2)->GetPositionVector());

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
