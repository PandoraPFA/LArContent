/**
 *  @file   larpandoracontent/LArObjects/LArTrackPfo.cc
 *
 *  @brief  Implementation of the lar track pfo class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArTrackPfo.h"

#include "Objects/CaloHit.h"

using namespace pandora;

namespace lar_content
{

LArTrackPfo::LArTrackPfo(const LArTrackPfoParameters &parameters) :
    ParticleFlowObject(parameters),
    m_trackStateVector(parameters.m_trackStateVector)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArTrackPfo::GetVertexPosition() const
{
    if (m_trackStateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return (*m_trackStateVector.begin()).GetPosition();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArTrackPfo::GetEndPosition() const
{
    if (m_trackStateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return (*m_trackStateVector.rbegin()).GetPosition();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArTrackPfo::GetVertexDirection() const
{
    if (m_trackStateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return (*m_trackStateVector.begin()).GetDirection();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArTrackPfo::GetEndDirection() const
{
    if (m_trackStateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return (*m_trackStateVector.rbegin()).GetDirection();
}

} // namespace lar_content
