/**
 *  @file   larpandoracontent/LArHelpers/LArCaloHitHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArCaloHitHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

float LArCaloHitHelper::GetClosestDistance(const CartesianVector &position, const CartesianPointVector &positionVector)
{
    return (position - LArCaloHitHelper::GetClosestPosition(position, positionVector)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArCaloHitHelper::GetClosestPosition(const CartesianVector &position, const CartesianPointVector &positionVector)
{
    bool found(false);
    CartesianVector closestPosition(0.f, 0.f, 0.f);
    float minDistSq(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : positionVector)
    {
        const float distSq((testPosition - position).GetMagnitudeSquared());

        if (distSq < minDistSq)
        {
            found = true;
            minDistSq = distSq;
            closestPosition = position;
        }
    }

    if (!found)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return closestPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArCaloHitHelper::IsInSameVolume(const CaloHit *const pPrevHit, const CaloHit *const pCurrentHit)
{
    if (!pPrevHit || !pCurrentHit)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LArCaloHit *const pPrevLArHit(dynamic_cast<const LArCaloHit *>(pPrevHit));
    const LArCaloHit *const pCurrentLArHit(dynamic_cast<const LArCaloHit *>(pCurrentHit));

    if (!pPrevLArHit || !pCurrentLArHit)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    unsigned int prevTPCID(pPrevLArHit->GetLArTPCVolumeId());
    unsigned int currentTPCID(pCurrentLArHit->GetLArTPCVolumeId());

    if (prevTPCID != currentTPCID)
        return false;

    unsigned int prevChildVolID(pPrevLArHit->GetDaughterVolumeId());
    unsigned int currentChildVolID(pCurrentLArHit->GetDaughterVolumeId());

    if (prevChildVolID != currentChildVolID)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
