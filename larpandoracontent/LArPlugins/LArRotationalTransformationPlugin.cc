/**
 *  @file   LArContent/LArPlugins/LArRotationalTransformationPlugin.cc
 *
 *  @brief  Implementation of the rotational transformation plugin class.
 *
 *  $Log: $
 */

#include "Objects/CartesianVector.h"

#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include <cmath>

namespace lar_content
{

using namespace pandora;

LArRotationalTransformationPlugin::LArRotationalTransformationPlugin(const double thetaU, const double thetaV, const double sigmaUVW) : 
    m_thetaU(thetaU),
    m_thetaV(thetaV),
    m_sigmaUVW(sigmaUVW),
    m_sinUminusV(std::sin(m_thetaU - m_thetaV)),
    m_sinUplusV(std::sin(m_thetaU + m_thetaV)),
    m_sinU(std::sin(m_thetaU)),
    m_sinV(std::sin(m_thetaV)),
    m_cosU(std::cos(m_thetaU)),
    m_cosV(std::cos(m_thetaV))
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

LArRotationalTransformationPlugin::~LArRotationalTransformationPlugin()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UVtoW(const double u, const double v) const
{
    return this->UVtoZ(u, v);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::VWtoU(const double v, const double w) const
{
    return (w * m_sinUplusV - v * m_sinU) / m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::WUtoV(const double w, const double u) const
{
    return (w * m_sinUplusV - u * m_sinV) / m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UVtoY(const double u, const double v) const
{
    return (v * m_cosU - u * m_cosV) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UVtoZ(const double u, const double v) const
{
    return (v * m_sinU + u * m_sinV) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::YZtoU(const double y, const double z) const
{
    return z * m_cosU - y * m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::YZtoV(const double y, const double z) const
{
    return z * m_cosV + y * m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::GetSigmaUVW() const
{
    return m_sigmaUVW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRotationalTransformationPlugin::GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV,
    const double sigmaW, double &y, double &z, double &chiSquared) const
{
    // Obtain expression for chi2, differentiate wrt y and z, set both results to zero and solve simultaneously.
    y = ((sigmaU * v * m_sinV) - (sigmaU * w * m_cosV * m_sinV) - (sigmaV * u * m_sinU) + (sigmaV * w * m_cosU * m_sinU) -
         (sigmaW * u * m_cosV * m_sinUplusV) + (sigmaW * v * m_cosU * m_sinUplusV)) /
        ((sigmaV * m_sinU * m_sinU) + (sigmaW * m_cosV * m_cosV * m_sinU * m_sinU) + (2. * sigmaW * m_cosU * m_cosV * m_sinU * m_sinV) +
         (sigmaU * m_sinV * m_sinV) + (sigmaW * m_cosU * m_cosU * m_sinV * m_sinV));

    z = ((sigmaU * w * m_sinV * m_sinV) + (sigmaV * w * m_sinU * m_sinU) +
         (sigmaW * u * m_sinV * m_sinUplusV) + (sigmaW * v * m_sinU * m_sinUplusV)) /
        ((sigmaV * m_sinU * m_sinU) + (sigmaW * m_cosV * m_cosV * m_sinU * m_sinU) + (2. * sigmaW * m_cosU * m_cosV * m_sinU * m_sinV) +
         (sigmaU * m_sinV * m_sinV) + (sigmaW * m_cosU * m_cosU * m_sinV * m_sinV));

    const double deltaU(u - LArRotationalTransformationPlugin::YZtoU(y, z));
    const double deltaV(v - LArRotationalTransformationPlugin::YZtoV(y, z));
    const double deltaW(w - z);
    chiSquared = ((deltaU * deltaU) / (sigmaU * sigmaU)) + ((deltaV * deltaV) / (sigmaV * sigmaV)) + ((deltaW * deltaW) / (sigmaW * sigmaW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRotationalTransformationPlugin::GetProjectedYZ(const PositionAndType &hitPositionAndType, const PositionAndType &fitPositionAndType1,
    const PositionAndType &fitPositionAndType2, const double sigmaHit, const double sigmaFit, double &y, double &z, double &chiSquared) const
{
    const HitType hitType(hitPositionAndType.second), fitType1(fitPositionAndType1.second), fitType2(fitPositionAndType2.second);

    if (((hitType == fitType1) || (hitType == fitType2) || (fitType1 == fitType2)) ||
        ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType)) ||
        ((TPC_VIEW_U != fitType1) && (TPC_VIEW_V != fitType1) && (TPC_VIEW_W != fitType1)) ||
        ((TPC_VIEW_U != fitType2) && (TPC_VIEW_V != fitType2) && (TPC_VIEW_W != fitType2)))
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const double u((TPC_VIEW_U == hitType) ? hitPositionAndType.first : (TPC_VIEW_U == fitType1) ? fitPositionAndType1.first : fitPositionAndType2.first);
    const double v((TPC_VIEW_V == hitType) ? hitPositionAndType.first : (TPC_VIEW_V == fitType1) ? fitPositionAndType1.first : fitPositionAndType2.first);
    const double w((TPC_VIEW_W == hitType) ? hitPositionAndType.first : (TPC_VIEW_W == fitType1) ? fitPositionAndType1.first : fitPositionAndType2.first);

    const double uPrime((TPC_VIEW_U == hitType) ? LArRotationalTransformationPlugin::VWtoU(v, w) : u);
    const double vPrime((TPC_VIEW_V == hitType) ? LArRotationalTransformationPlugin::WUtoV(w, u) : v);
    const double yInput(LArRotationalTransformationPlugin::UVtoY(uPrime, vPrime));
    const double zInput(LArRotationalTransformationPlugin::UVtoZ(uPrime, vPrime));

    const double position((TPC_VIEW_U == hitType) ? u : (TPC_VIEW_V == hitType) ? v : w);
    const double fitPosition((TPC_VIEW_U == hitType) ? LArRotationalTransformationPlugin::VWtoU(v, w) : (TPC_VIEW_V == hitType) ? LArRotationalTransformationPlugin::WUtoV(w, u) : LArRotationalTransformationPlugin::UVtoW(u, v));
    const double unitY((TPC_VIEW_U == hitType) ? -m_cosU : (TPC_VIEW_V == hitType) ? m_cosV : 0.);
    const double unitZ((TPC_VIEW_U == hitType) ? m_sinU : (TPC_VIEW_V == hitType) ? m_sinV : 1.);  

    y = yInput + unitY * (position - fitPosition);
    z = zInput + unitZ * (position - fitPosition);

    const double sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
    const double sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
    const double sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);
    const double deltaU(u - LArRotationalTransformationPlugin::YZtoU(y, z));
    const double deltaV(v - LArRotationalTransformationPlugin::YZtoV(y, z));
    const double deltaW(w - z);
    chiSquared = ((deltaU * deltaU) / (sigmaU * sigmaU)) + ((deltaV * deltaV) / (sigmaV * sigmaV)) + ((deltaW * deltaW) / (sigmaW * sigmaW));
}

} // namespace lar_content
