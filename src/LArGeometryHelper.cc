/**
 *  @file   LArGeometryHelper.cc
 * 
 *  @brief  Implementation of the geometry helper class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "LArGeometryHelper.h"

using namespace pandora;

namespace lar
{

float LArGeometryHelper::MergeTwoPositions(const HitType view1, const HitType view2, const float position1, const float position2)
{
    float U(0.f), V(0.f), W(0.f);
    HitType newView(INNER_DETECTOR);

    if (view1 == VIEW_U && view2 == VIEW_V)
    {
        U = position1; V = position2; newView = VIEW_W;
    }

    if (view1 == VIEW_V && view2 == VIEW_U)
    {
        V = position1; U = position2; newView = VIEW_W;
    }

    if (view1 == VIEW_W && view2 == VIEW_U)
    {
        W = position1; U = position2; newView = VIEW_V;
    }

    if (view1 == VIEW_U && view2 == VIEW_W)
    {
        U = position1; W = position2; newView = VIEW_V;
    }

    if (view1 == VIEW_V && view2 == VIEW_W)
    {
        V = position1; W = position2; newView = VIEW_U;
    }

    if (view1 == VIEW_W && view2 == VIEW_V)
    {
        W = position1; V = position2; newView = VIEW_U;
    }

    if (newView == VIEW_W)
        return (U * m_sinV + V * m_sinU - m_H * m_sinU * m_sinV) / m_sinUV;

    if (newView == VIEW_U)
        return (W * m_sinUV - V * m_sinU + m_H * m_sinU * m_sinV) / m_sinV;

    if (newView == VIEW_V)
        return (W * m_sinUV - U * m_sinV + m_H * m_sinU * m_sinV) / m_sinU;

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::MergeTwoDirections(const HitType view1, const HitType view2, const CartesianVector &direction1,
    const CartesianVector &direction2)
{
    // x components must be consistent
    if (direction1.GetX() * direction2.GetX() < 0.f)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // scale x components to a common value
    const float s1((std::fabs(direction1.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction2.GetX()) : 1.f);
    const float s2((std::fabs(direction2.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction1.GetX()) : 1.f);

    float pX(s1 * direction1.GetX()), pU(0.f), pV(0.f), pW(0.f);
    HitType newView(INNER_DETECTOR);

    if (view1 == VIEW_U && view2 == VIEW_V)
    {
        pU = s1 * direction1.GetZ(); pV = s2 * direction2.GetZ(); newView = VIEW_W;
    }

    if (view1 == VIEW_V && view2 == VIEW_U)
    {
        pV = s1 * direction1.GetZ(); pU = s2 * direction2.GetZ(); newView = VIEW_W;
    }

    if (view1 == VIEW_W && view2 == VIEW_U)
    {
        pW = s1 * direction1.GetZ(); pU = s2 * direction2.GetZ(); newView = VIEW_V;
    }

    if (view1 == VIEW_U && view2 == VIEW_W)
    {
        pU = s1 * direction1.GetZ(); pW = s2 * direction2.GetZ(); newView = VIEW_V;
    }

    if (view1 == VIEW_V && view2 == VIEW_W)
    {
        pV = s1 * direction1.GetZ(); pW = s2 * direction2.GetZ(); newView = VIEW_U;
    }

    if (view1 == VIEW_W && view2 == VIEW_V)
    {
        pW = s1 * direction1.GetZ(); pV = s2 * direction2.GetZ(); newView = VIEW_U;//
    }

    if (newView == VIEW_W)
        return CartesianVector(pX, 0.f, (pU * m_sinV + pV * m_sinU) / m_sinUV).GetUnitVector();

    if (newView == VIEW_U)
        return CartesianVector(pX, 0.f, (pW * m_sinUV - pV * m_sinU) / m_sinV).GetUnitVector();

    if (newView == VIEW_V)
        return CartesianVector(pX, 0.f, (pW * m_sinUV - pU * m_sinV) / m_sinU).GetUnitVector();

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::m_thetaU = M_PI / 3.f;
float LArGeometryHelper::m_thetaV = M_PI / 3.f;
float LArGeometryHelper::m_H = 233.f;

float LArGeometryHelper::m_sinU  = std::sin(LArGeometryHelper::m_thetaU);
float LArGeometryHelper::m_sinV  = std::sin(LArGeometryHelper::m_thetaV);
float LArGeometryHelper::m_sinUV = std::sin(LArGeometryHelper::m_thetaU + LArGeometryHelper::m_thetaV);

StatusCode LArGeometryHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    float thetaU_degrees = 60.f;
    if (STATUS_CODE_SUCCESS == XmlHelper::ReadValue(xmlHandle, "ThetaU", thetaU_degrees))
    {
        m_thetaU = thetaU_degrees * (M_PI / 180.f);
        m_sinU = std::sin(m_thetaU);
    }

    float thetaV_degrees = 60.f;
    if (STATUS_CODE_SUCCESS == XmlHelper::ReadValue(xmlHandle, "ThetaU", thetaV_degrees))
    {
        m_thetaV = thetaV_degrees * (M_PI / 180.f);
        m_sinV = std::sin(m_thetaV);
    }

    m_sinUV = std::sin(m_thetaU + m_thetaV);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Height", m_H));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
