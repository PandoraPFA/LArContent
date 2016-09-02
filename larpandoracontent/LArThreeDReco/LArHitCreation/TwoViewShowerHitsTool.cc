/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.cc
 *
 *  @brief  Implementation of the two view shower hits tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h"

using namespace pandora;

namespace lar_content
{

void TwoViewShowerHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const CaloHitList &caloHitList1, const CaloHitList &caloHitList2,
    CartesianVector &position3D, float &chiSquared) const
{
    if (!caloHitList1.empty() && caloHitList2.empty())
    {
        this->GetThreeDPosition(pCaloHit2D, caloHitList1, position3D, chiSquared);
    }
    else if (caloHitList1.empty() && !caloHitList2.empty())
    {
        this->GetThreeDPosition(pCaloHit2D, caloHitList2, position3D, chiSquared);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewShowerHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const CaloHitList &caloHitList, CartesianVector &position3D,
    float &chiSquared) const
{
    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType hitType((*caloHitList.begin())->GetHitType());

    if (pCaloHit2D->GetHitType() == hitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    double Sqz(0.), Sqx(0.), Sq(0.);

    for (CaloHitList::const_iterator iter = caloHitList.begin(), iterEnd = caloHitList.end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pCaloHit = *iter;

        Sqx += pCaloHit->GetMipEquivalentEnergy() * pCaloHit->GetPositionVector().GetX();
        Sqz += pCaloHit->GetMipEquivalentEnergy() * pCaloHit->GetPositionVector().GetZ();
        Sq  += pCaloHit->GetMipEquivalentEnergy();
    }

    if (Sq < std::numeric_limits<double>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector position(static_cast<float>(Sqx / Sq), 0.f, static_cast<float>(Sqz / Sq));

    this->GetPosition3D(pCaloHit2D, hitType, position, position3D, chiSquared);
}

} // namespace lar_content
