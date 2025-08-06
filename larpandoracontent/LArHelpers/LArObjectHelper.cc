/**
 *  @file   larpandoracontent/LArHelpers/LArObjectHelper.cc
 *
 *  @brief  Implementation of the object helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArObjectHelper.h"

using namespace pandora;

namespace lar_content
{

template <>
CartesianVector LArObjectHelper::TypeAdaptor::GetPosition(const CartesianVector &t)
{
    return t;
}

template <>
CartesianVector LArObjectHelper::TypeAdaptor::GetPosition(const CaloHit *const &pT)
{
    return pT->GetPositionVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
const CaloHit *LArObjectHelper::TypeAdaptor::GetCaloHit(const CartesianVector &, [[maybe_unused]] const bool retSelf)
{
    return nullptr;
}

template <>
const CaloHit *LArObjectHelper::TypeAdaptor::GetCaloHit(const CaloHit *const &pCaloHit3D, const bool retSelf)
{
    if (retSelf)
        return static_cast<const CaloHit *>(pCaloHit3D);
    else
        return static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress());
}

} // namespace lar_content
