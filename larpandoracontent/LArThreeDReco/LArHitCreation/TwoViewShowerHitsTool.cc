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

void TwoViewShowerHitsTool::GetShowerHit3D(const CaloHitVector &caloHitVector1, const CaloHitVector &caloHitVector2, ProtoHit &protoHit) const
{
    if (!caloHitVector1.empty() && caloHitVector2.empty())
    {
        this->GetShowerHit3D(caloHitVector1, protoHit);
    }
    else if (caloHitVector1.empty() && !caloHitVector2.empty())
    {
        this->GetShowerHit3D(caloHitVector2, protoHit);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewShowerHitsTool::GetShowerHit3D(const CaloHitVector &caloHitVector, ProtoHit &protoHit) const
{
    if (caloHitVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CaloHit *const pCaloHit2D(protoHit.GetParentCaloHit2D());
    const HitType hitType(caloHitVector.at(0)->GetHitType());

    if (pCaloHit2D->GetHitType() == hitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    double Sqz(0.), Sqx(0.), Sq(0.);

    for (const CaloHit *const pCaloHit : caloHitVector)
    {
        Sqx += pCaloHit->GetMipEquivalentEnergy() * pCaloHit->GetPositionVector().GetX();
        Sqz += pCaloHit->GetMipEquivalentEnergy() * pCaloHit->GetPositionVector().GetZ();
        Sq  += pCaloHit->GetMipEquivalentEnergy();
    }

    if (Sq < std::numeric_limits<double>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector position(static_cast<float>(Sqx / Sq), 0.f, static_cast<float>(Sqz / Sq));
    this->GetBestPosition3D(hitType, position, protoHit);

    // ATTN For this treatment of showers, deltaX term in chi2 is important and is added here
    const double deltaX(pCaloHit2D->GetPositionVector().GetX() - position.GetX());
    const double chi2X((deltaX * deltaX) / this->GetSigmaX2());
    protoHit.SetPosition3D(protoHit.GetPosition3D(), protoHit.GetChi2() + chi2X);
}

} // namespace lar_content
