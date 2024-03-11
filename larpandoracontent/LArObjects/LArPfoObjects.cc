/**
 *  @file   larpandoracontent/LArObjects/LArPfoObjects.cc
 *
 *  @brief  Implementation of lar pfo objects.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include "Objects/CaloHit.h"

using namespace pandora;

namespace lar_content
{

LArTrackState::LArTrackState(const CartesianVector &position, const CartesianVector &direction, const CaloHit *const pCaloHit) :
    TrackState(position, direction),
    m_pCaloHit(pCaloHit)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArTrackState::LArTrackState(const CartesianVector &position, const CartesianVector &direction) :
    TrackState(position, direction),
    m_pCaloHit(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArTrackState::GetDirection() const
{
    return this->GetMomentum();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *LArTrackState::GetCaloHit() const
{
    if (m_pCaloHit)
        return m_pCaloHit;

    throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArShowerPCA::LArShowerPCA(const CartesianVector &centroid, const CartesianVector &primaryAxis, const CartesianVector &secondaryAxis,
    const CartesianVector &tertiaryAxis, const CartesianVector &eigenValues) :
    m_centroid(centroid),
    m_primaryAxis(primaryAxis),
    m_secondaryAxis(secondaryAxis),
    m_tertiaryAxis(tertiaryAxis),
    m_eigenValues(eigenValues),
    m_axisLengths(((eigenValues.GetX() > std::numeric_limits<float>::epsilon()) ? 6.f * std::sqrt(eigenValues.GetX()) : 0.f),
        ((eigenValues.GetY() > std::numeric_limits<float>::epsilon()) ? 6.f * std::sqrt(eigenValues.GetY()) : 0.f),
        ((eigenValues.GetZ() > std::numeric_limits<float>::epsilon()) ? 6.f * std::sqrt(eigenValues.GetZ()) : 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArShowerPCA::GetCentroid() const
{
    return m_centroid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArShowerPCA::GetPrimaryAxis() const
{
    return m_primaryAxis;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArShowerPCA::GetSecondaryAxis() const
{
    return m_secondaryAxis;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArShowerPCA::GetTertiaryAxis() const
{
    return m_tertiaryAxis;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArShowerPCA::GetEigenValues() const
{
    return m_eigenValues;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArShowerPCA::GetAxisLengths() const
{
    return m_axisLengths;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArShowerPCA::GetPrimaryLength() const
{
    return m_axisLengths.GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArShowerPCA::GetSecondaryLength() const
{
    return m_axisLengths.GetY();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArShowerPCA::GetTertiaryLength() const
{
    return m_axisLengths.GetZ();
}

} // namespace lar_content
