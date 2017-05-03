/**
 *  @file   larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.cc
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
    const double sigmaU2(sigmaU * sigmaU), sigmaV2(sigmaV * sigmaV), sigmaW2(sigmaW * sigmaW);

    // Obtain expression for chi2, differentiate wrt y and z, set both results to zero and solve simultaneously. Here just paste-in result.
    y = ((sigmaU2 * v * m_sinV) - (sigmaU2 * w * m_cosV * m_sinV) - (sigmaV2 * u * m_sinU) + (sigmaV2 * w * m_cosU * m_sinU) -
         (sigmaW2 * u * m_cosV * m_sinUplusV) + (sigmaW2 * v * m_cosU * m_sinUplusV)) /
        ((sigmaV2 * m_sinU * m_sinU) + (sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU) + (2. * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV) +
         (sigmaU2 * m_sinV * m_sinV) + (sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV));

    z = ((sigmaU2 * w * m_sinV * m_sinV) + (sigmaV2 * w * m_sinU * m_sinU) +
         (sigmaW2 * u * m_sinV * m_sinUplusV) + (sigmaW2 * v * m_sinU * m_sinUplusV)) /
        ((sigmaV2 * m_sinU * m_sinU) + (sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU) + (2. * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV) +
         (sigmaU2 * m_sinV * m_sinV) + (sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV));

    const double deltaU(u - LArRotationalTransformationPlugin::YZtoU(y, z));
    const double deltaV(v - LArRotationalTransformationPlugin::YZtoV(y, z));
    const double deltaW(w - z);
    chiSquared = ((deltaU * deltaU) / sigmaU2) + ((deltaV * deltaV) / sigmaV2) + ((deltaW * deltaW) / sigmaW2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRotationalTransformationPlugin::GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV,
    const double sigmaW, const double uFit, const double vFit, const double wFit, const double sigmaFit, double &y, double &z, double &chiSquared) const
{
    const double sigmaU2(sigmaU * sigmaU), sigmaV2(sigmaV * sigmaV), sigmaW2(sigmaW * sigmaW), sigmaFit2(sigmaFit * sigmaFit);

    // Obtain expression for chi2, differentiate wrt y and z, set both results to zero and solve simultaneously. Here just paste-in result.
    y = (-uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_sinU - uFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_sinU - sigmaV2 * sigmaW2 * sigmaFit2 * u * m_sinU -
        sigmaV2 * sigmaFit2 * sigmaFit2 * u * m_sinU + wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinU +
        wFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_sinU + sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosU * m_sinU +
        sigmaV2 * sigmaFit2 * sigmaFit2 * w * m_cosU * m_sinU + vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinU +
        vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU +
        sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosU * m_cosV * m_sinU +
        sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosU * m_cosV * m_sinU - uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinU -
        uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU - sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosV * m_cosV * m_sinU -
        sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosV * m_cosV * m_sinU + vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_sinV +
        vFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * v * m_sinV + sigmaU2 * sigmaFit2 * sigmaFit2 * v * m_sinV +
        vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinV + vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV +
        sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosU * m_cosU * m_sinV + sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosU * m_cosU * m_sinV -
        wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinV - wFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_sinV -
        sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosV * m_sinV - sigmaU2 * sigmaFit2 * sigmaFit2 * w * m_cosV * m_sinV -
        uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinV -
        uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinV -
        sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosU * m_cosV * m_sinV -
        sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosU * m_cosV * m_sinV)
        /
        (sigmaU2 * sigmaV2 * sigmaW2 * m_sinU * m_sinU +
        sigmaU2 * sigmaV2 * sigmaFit2 * m_sinU * m_sinU + sigmaV2 * sigmaW2 * sigmaFit2 * m_sinU * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * m_sinU * m_sinU +
        sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
        sigmaV2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
        2 * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV +
        2 * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV +
        2 * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV +
        2 * sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV + sigmaU2 * sigmaV2 * sigmaW2 * m_sinV * m_sinV +
        sigmaU2 * sigmaV2 * sigmaFit2 * m_sinV * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * m_sinV * m_sinV + sigmaU2 * sigmaFit2 * sigmaFit2 * m_sinV * m_sinV +
        sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV +
        sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV);

    z = (wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_sinU * m_sinU + wFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_sinU * m_sinU +
        sigmaU2 * sigmaV2 * sigmaFit2 * w * m_sinU * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * w * m_sinU * m_sinU +
        vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinU * m_sinU + vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosV * m_sinU * m_sinU +
        sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosV * m_sinU * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosV * m_sinU * m_sinU +
        vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinU * m_sinV +
        vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_sinU * m_sinV +
        sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosU * m_sinU * m_sinV +
        sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosU * m_sinU * m_sinV +
        uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinU * m_sinV +
        uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_sinU * m_sinV +
        sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosV * m_sinU * m_sinV +
        sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosV * m_sinU * m_sinV + wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_sinV * m_sinV +
        wFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_sinV * m_sinV + sigmaU2 * sigmaV2 * sigmaFit2 * w * m_sinV * m_sinV +
        sigmaU2 * sigmaFit2 * sigmaFit2 * w * m_sinV * m_sinV + uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinV * m_sinV +
        uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_sinV * m_sinV + sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosU * m_sinV * m_sinV +
        sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosU * m_sinV * m_sinV)
        /
        (sigmaU2 * sigmaV2 * sigmaW2 * m_sinU * m_sinU +
        sigmaU2 * sigmaV2 * sigmaFit2 * m_sinU * m_sinU + sigmaV2 * sigmaW2 * sigmaFit2 * m_sinU * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * m_sinU * m_sinU +
        sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
        sigmaV2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
        2 * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV +
        2 * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV +
        2 * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV +
        2 * sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV + sigmaU2 * sigmaV2 * sigmaW2 * m_sinV * m_sinV +
        sigmaU2 * sigmaV2 * sigmaFit2 * m_sinV * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * m_sinV * m_sinV + sigmaU2 * sigmaFit2 * sigmaFit2 * m_sinV * m_sinV +
        sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV +
        sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV);

    const double outputU(LArRotationalTransformationPlugin::YZtoU(y, z));
    const double outputV(LArRotationalTransformationPlugin::YZtoV(y, z));

    const double deltaU(u - outputU), deltaV(v - outputV), deltaW(w - z);
    const double deltaUFit(uFit - outputU), deltaVFit(vFit - outputV), deltaWFit(wFit - z);

    chiSquared = ((deltaU * deltaU) / sigmaU2) + ((deltaV * deltaV) / sigmaV2) + ((deltaW * deltaW) / sigmaW2) +
        ((deltaUFit * deltaUFit) / sigmaFit2) + ((deltaVFit * deltaVFit) / sigmaFit2) + ((deltaWFit * deltaWFit) / sigmaFit2);
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
