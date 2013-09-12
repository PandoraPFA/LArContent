/**
 *  @file   LArContent/src/LArHelpers/LArGeometryHelper.cc
 * 
 *  @brief  Implementation of the geometry helper class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar
{

float LArGeometryHelper::MergeTwoPositions(const HitType view1, const HitType view2, const float position1, const float position2)
{
    if (view1 == view2)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((view1 == VIEW_U) && (view2 == VIEW_V))
    {
        return LArGeometryHelper::UVtoW(position1, position2);
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_U))
    {
        return LArGeometryHelper::VUtoW(position1, position2);
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_U))
    {
        return LArGeometryHelper::WUtoV(position1, position2);
    }

    if ((view1 == VIEW_U) && (view2 == VIEW_W))
    {
        return LArGeometryHelper::UWtoV(position1, position2);
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_W))
    {
        return LArGeometryHelper::VWtoU(position1, position2);
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_V))
    {
        return LArGeometryHelper::WVtoU(position1, position2);
    }

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::MergeTwoDirections(const HitType view1, const HitType view2, const CartesianVector &direction1,
    const CartesianVector &direction2)
{
    // x components must be consistent
    if (direction1.GetX() * direction2.GetX() < 0.f)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // TODO: Can do better than this...

    // scale x components to a common value
    const float s1((std::fabs(direction1.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction2.GetX()) : 1.f);
    const float s2((std::fabs(direction2.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction1.GetX()) : 1.f);

    float pX(s1 * direction1.GetX()), pU(0.f), pV(0.f), pW(0.f);
    HitType newView(INNER_DETECTOR);

    if ((view1 == VIEW_U) && (view2 == VIEW_V))
    {
        pU = s1 * direction1.GetZ(); pV = s2 * direction2.GetZ(); newView = VIEW_W;
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_U))
    {
        pV = s1 * direction1.GetZ(); pU = s2 * direction2.GetZ(); newView = VIEW_W;
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_U))
    {
        pW = s1 * direction1.GetZ(); pU = s2 * direction2.GetZ(); newView = VIEW_V;
    }

    if ((view1 == VIEW_U) && (view2 == VIEW_W))
    {
        pU = s1 * direction1.GetZ(); pW = s2 * direction2.GetZ(); newView = VIEW_V;
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_W))
    {
        pV = s1 * direction1.GetZ(); pW = s2 * direction2.GetZ(); newView = VIEW_U;
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_V))
    {
        pW = s1 * direction1.GetZ(); pV = s2 * direction2.GetZ(); newView = VIEW_U;
    }

    if (newView == VIEW_W)
        return CartesianVector(pX, 0.f, LArGeometryHelper::PUPVtoPW(pU, pV)).GetUnitVector();
    
    if (newView == VIEW_U)
        return CartesianVector(pX, 0.f, LArGeometryHelper::PWPVtoPU(pW, pV)).GetUnitVector();
    
    if (newView == VIEW_V)
        return CartesianVector(pX, 0.f, LArGeometryHelper::PWPUtoPV(pW, pU)).GetUnitVector();

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions(const HitType view1, const HitType view2, const CartesianVector &position1, const CartesianVector &position2, CartesianVector &position3, float& chiSquared)
{
    if (view1 == view2)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float X3((position1.GetX() + position2.GetX() ) / 2.f);
    const float Z1(position1.GetZ());
    const float Z2(position2.GetZ());
    const float Z3(LArGeometryHelper::MergeTwoPositions(view1, view2, Z1, Z2));

    position3.SetValues(X3, 0.f, Z3);
    chiSquared = ((X3 - position1.GetX()) * (X3 - position1.GetX()) + (X3 - position2.GetX()) * (X3 - position2.GetX())) / (m_sigmaUVW * m_sigmaUVW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions(const HitType view1, const HitType view2, const CartesianVector &position1, const CartesianVector &position2, CartesianVector &outputU, CartesianVector &outputV, CartesianVector &outputW, float& chiSquared)
{
    if (view1 == view2)
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector position3(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(view1, view2, position1, position2, position3, chiSquared);

    float aveU(0.f), aveV(0.f), aveW(0.f);
    const float Z1(position1.GetZ()), Z2(position2.GetZ()), Z3(position3.GetZ()), aveX(position3.GetX());

    if ((view1 == VIEW_U) && (view2 == VIEW_V))
    {
        aveU = Z1; aveV = Z2; aveW = Z3;
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_W))
    {
        aveV = Z1; aveW = Z2; aveU = Z3;
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_U))
    {
        aveW = Z1; aveU = Z2; aveV = Z3;
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_U))
    {
        aveV = Z1; aveU = Z2; aveW = Z3;
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_V))
    {
        aveW = Z1; aveV = Z2; aveU = Z3;
    }

    if ((view1 == VIEW_U) && (view2 == VIEW_W))
    {
        aveU = Z1; aveW = Z2; aveV = Z3;
    }

    outputU.SetValues( aveX, 0.f, aveU );
    outputV.SetValues( aveX, 0.f, aveV );
    outputW.SetValues( aveX, 0.f, aveW );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions(const HitType view1, const HitType view2, const HitType view3, const CartesianVector &position1,
    const CartesianVector &position2, const CartesianVector &position3, CartesianVector &outputU, CartesianVector &outputV,
    CartesianVector &outputW, float &chiSquared)
{
    if ((view1 == view2) || (view2 == view3) || (view3 == view1))
        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((view1 == VIEW_U) && (view2 == VIEW_V))
    {
        return LArGeometryHelper::MergeThreePositions(position1, position2, position3, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_W))
    {
        return LArGeometryHelper::MergeThreePositions(position3, position1, position2, outputU, outputV, outputW, chiSquared );
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_U))
    {
        return LArGeometryHelper::MergeThreePositions(position2, position3, position1, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == VIEW_V) && (view2 == VIEW_U))
    {
        return LArGeometryHelper::MergeThreePositions(position2, position1, position3, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == VIEW_W) && (view2 == VIEW_V))
    {
        return LArGeometryHelper::MergeThreePositions(position3, position2, position1, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == VIEW_U) && (view2 == VIEW_W))
    {
        return LArGeometryHelper::MergeThreePositions(position1, position3, position2, outputU, outputV, outputW, chiSquared);
    }

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions(const CartesianVector &positionU, const CartesianVector &positionV, const CartesianVector &positionW,
    CartesianVector &outputU, CartesianVector &outputV, CartesianVector &outputW, float &chiSquared)
{
    const float YfromUV(LArGeometryHelper::UVtoY(positionU.GetZ(), positionV.GetZ()));
    const float ZfromUV(LArGeometryHelper::UVtoZ(positionU.GetZ(), positionV.GetZ()));

    const float aveX((positionU.GetX() + positionV.GetX() + positionW.GetX()) / 3.f);
    const float aveY(YfromUV);
    const float aveW((positionW.GetZ() + 2.f * ZfromUV ) / 3.f);

    const float aveU(LArGeometryHelper::YZtoU(aveY, aveW));
    const float aveV(LArGeometryHelper::YZtoV(aveY, aveW));

    outputU.SetValues(aveX, 0.f, aveU);
    outputV.SetValues(aveX, 0.f, aveV);
    outputW.SetValues(aveX, 0.f, aveW);

    chiSquared = ((outputU.GetX() - positionU.GetX()) * (outputU.GetX() - positionU.GetX()) +
                  (outputV.GetX() - positionV.GetX()) * (outputV.GetX() - positionV.GetX()) +
                  (outputW.GetX() - positionW.GetX()) * (outputW.GetX() - positionW.GetX()) +
                  (outputU.GetZ() - positionU.GetZ()) * (outputU.GetZ() - positionU.GetZ()) +
                  (outputV.GetZ() - positionV.GetZ()) * (outputV.GetZ() - positionV.GetZ()) +
                  (outputW.GetZ() - positionW.GetZ()) * (outputW.GetZ() - positionW.GetZ())) / (m_sigmaUVW * m_sigmaUVW);
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions3D(const HitType view1, const HitType view2, const CartesianVector &position1,
    const CartesianVector &position2, CartesianVector &position3D, float &chiSquared)
{
    CartesianVector positionU(0.f, 0.f, 0.f), positionV(0.f, 0.f, 0.f), positionW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(view1, view2, position1, position2, positionU, positionV, positionW, chiSquared);

    position3D.SetValues(positionW.GetX(), LArGeometryHelper::UVtoY(positionU.GetZ(),positionV.GetZ()),
      LArGeometryHelper::UVtoZ(positionU.GetZ(),positionV.GetZ()) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions3D(const HitType view1, const HitType view2, const HitType view3, const CartesianVector &position1,
    const CartesianVector &position2, const CartesianVector &position3, CartesianVector &position3D, float &chiSquared)
{
    CartesianVector positionU(0.f, 0.f, 0.f), positionV(0.f, 0.f, 0.f), positionW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeThreePositions(view1, view2, view3, position1, position2, position3, positionU, positionV, positionW, chiSquared);

    position3D.SetValues(positionW.GetX(), LArGeometryHelper::UVtoY(positionU.GetZ(),positionV.GetZ()),
        LArGeometryHelper::UVtoZ(positionU.GetZ(),positionV.GetZ()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::ProjectPosition(const CartesianVector &position3D, const HitType view)
{
    if (view == VIEW_U)
    {
        return CartesianVector(position3D.GetX(), 0.0, LArGeometryHelper::YZtoU(position3D.GetY(), position3D.GetZ()));
    }

    else if (view == VIEW_V)
    {
        return CartesianVector(position3D.GetX(), 0.0, LArGeometryHelper::YZtoV(position3D.GetY(), position3D.GetZ()));
    }

    else if (view == VIEW_W)
    {
        return CartesianVector(position3D.GetX(), 0.0, position3D.GetZ());
    }

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
CartesianVector LArGeometryHelper::ProjectDirection(const CartesianVector &direction3D, const HitType view)
{
    if (view == VIEW_U)
    {
        return CartesianVector(direction3D.GetX(), 0.f, LArGeometryHelper::PYPZtoPU(direction3D.GetY(), direction3D.GetZ())).GetUnitVector();
    }

    else if (view == VIEW_V)
    {
        return CartesianVector(direction3D.GetX(), 0.f, LArGeometryHelper::PYPZtoPV(direction3D.GetY(), direction3D.GetZ())).GetUnitVector();
    }

    else if (view == VIEW_W)
    {
        return CartesianVector(direction3D.GetX(), 0.f, direction3D.GetZ()).GetUnitVector();
    }

    throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::VUtoW(const float v, const float u)
{
    return LArGeometryHelper::UVtoW(u, v);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::WVtoU(const float w, const float v)
{
    return LArGeometryHelper::VWtoU(v, w);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::UWtoV(const float u, const float w)
{
    return LArGeometryHelper::WUtoV(w, u);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::UVtoW(const float u, const float v)
{
    return LArGeometryHelper::UVtoZ(u, v);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::VWtoU(const float v, const float w)
{
    return (w * m_sinUplusV - v * m_sinU + m_H * m_sinU * m_sinV) / m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::WUtoV(const float w, const float u)
{
    return (w * m_sinUplusV - u * m_sinV + m_H * m_sinU * m_sinV) / m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::UVtoY(const float u, const float v)
{
    return (u * m_cosV - v * m_cosU - 0.5f * m_H * m_sinUminusV) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::UVtoZ(const float u, const float v)
{
    return (u * m_sinV + v * m_sinU - m_H * m_sinU * m_sinV) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::YZtoU(const float y, const float z)
{
    return z * m_cosU + (y + 0.5 * m_H) * m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::YZtoV(const float y, const float z)
{
    return z * m_cosV - (y - 0.5 * m_H) * m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::PVPUtoPW(const float pv, const float pu)
{
    return LArGeometryHelper::PUPVtoPW(pu, pv);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::PWPVtoPU(const float pw, const float pv)
{
    return LArGeometryHelper::PVPWtoPU(pv, pw);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::PUPWtoPV(const float pu, const float pw)
{
    return LArGeometryHelper::PWPUtoPV(pw, pu);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::PUPVtoPW(const float pu, const float pv)
{
    return (pu * m_sinV + pv * m_sinU) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------
   
float LArGeometryHelper::PVPWtoPU(const float pv, const float pw)
{
    return (pw * m_sinUplusV - pv * m_sinU) / m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------
  
float LArGeometryHelper::PWPUtoPV(const float pw, const float pu)
{
    return (pw * m_sinUplusV - pu * m_sinV) / m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::PYPZtoPU(const float py, const float pz)
{
    return pz * m_cosU + py * m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::PYPZtoPV(const float py, const float pz)
{
    return pz * m_cosV - py * m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::m_thetaU = M_PI / 3.f;
float LArGeometryHelper::m_thetaV = M_PI / 3.f;
float LArGeometryHelper::m_H = 233.f;

float LArGeometryHelper::m_cosU  = std::cos(LArGeometryHelper::m_thetaU);
float LArGeometryHelper::m_cosV  = std::cos(LArGeometryHelper::m_thetaV);

float LArGeometryHelper::m_sinU  = std::sin(LArGeometryHelper::m_thetaU);
float LArGeometryHelper::m_sinV  = std::sin(LArGeometryHelper::m_thetaV);

float LArGeometryHelper::m_sinUplusV = std::sin(LArGeometryHelper::m_thetaU + LArGeometryHelper::m_thetaV);
float LArGeometryHelper::m_sinUminusV = std::sin(LArGeometryHelper::m_thetaU - LArGeometryHelper::m_thetaV);

float LArGeometryHelper::m_sigmaUVW = 1.f;

StatusCode LArGeometryHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    float thetaU_degrees = 60.f;
    if (STATUS_CODE_SUCCESS == XmlHelper::ReadValue(xmlHandle, "ThetaU", thetaU_degrees))
    {
        m_thetaU = thetaU_degrees * (M_PI / 180.f);
        m_sinU = std::sin(m_thetaU);
        m_cosU = std::cos(m_thetaU);
    }

    float thetaV_degrees = 60.f;
    if (STATUS_CODE_SUCCESS == XmlHelper::ReadValue(xmlHandle, "ThetaU", thetaV_degrees))
    {
        m_thetaV = thetaV_degrees * (M_PI / 180.f);
        m_sinV = std::sin(m_thetaV);
        m_cosV = std::cos(m_thetaV);
    }

    m_sinUplusV = std::sin(m_thetaU + m_thetaV);
    m_sinUminusV = std::sin(m_thetaU - m_thetaV);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Height", m_H));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Resolution", m_sigmaUVW));

    if (std::fabs(m_sigmaUVW) < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_INVALID_PARAMETER;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
