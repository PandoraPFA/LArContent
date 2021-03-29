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
const CaloHit *LArObjectHelper::TypeAdaptor::GetCaloHit(const CartesianVector &)
{
    return nullptr;
}

template <>
const CaloHit *LArObjectHelper::TypeAdaptor::GetCaloHit(const CaloHit *const &pCaloHit3D)
{
    const CaloHit *const pCaloHit2D = static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress());
    return pCaloHit2D;
}

} // namespace lar_content
