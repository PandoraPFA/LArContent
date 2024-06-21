/**
 *  @file   larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.cc
 *
 *  @brief  Implementation of the rotational transformation plugin class.
 *
 *  $Log: $
 */

#include "Geometry/LArTPC.h"

#include "Helpers/XmlHelper.h"

#include "Managers/GeometryManager.h"

#include "Objects/CartesianVector.h"

#include "Pandora/Pandora.h"

#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include <cmath>

namespace lar_content
{

using namespace pandora;

LArRotationalTransformationPlugin::LArRotationalTransformationPlugin() :
    m_thetaU(0.),
    m_thetaV(0.),
    m_thetaW(0.),
    m_sinU(0.),
    m_sinV(0.),
    m_sinW(0.),
    m_cosU(0.),
    m_cosV(0.),
    m_cosW(0.),
    m_sinVminusU(0.),
    m_sinWminusV(0.),
    m_sinUminusW(0.),
    m_maxAngularDiscrepancyU(0.03),
    m_maxAngularDiscrepancyV(0.03),
    m_maxAngularDiscrepancyW(0.03),
    m_maxSigmaDiscrepancy(0.01)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UVtoW(const double u, const double v) const
{
    return (-1. * (u * m_sinWminusV + v * m_sinUminusW) / m_sinVminusU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::VWtoU(const double v, const double w) const
{
    return (-1. * (v * m_sinUminusW + w * m_sinVminusU) / m_sinWminusV);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::WUtoV(const double w, const double u) const
{
    return (-1. * (u * m_sinWminusV + w * m_sinVminusU) / m_sinUminusW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UVtoY(const double u, const double v) const
{
    return ((u * m_cosV - v * m_cosU) / m_sinVminusU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UVtoZ(const double u, const double v) const
{
    return ((u * m_sinV - v * m_sinU) / m_sinVminusU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UWtoY(const double u, const double w) const
{
    return ((w * m_cosU - u * m_cosW) / m_sinUminusW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::UWtoZ(const double u, const double w) const
{
    return ((w * m_sinU - u * m_sinW) / m_sinUminusW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::VWtoY(const double v, const double w) const
{
    return ((v * m_cosW - w * m_cosV) / m_sinWminusV);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::VWtoZ(const double v, const double w) const
{
    return ((v * m_sinW - w * m_sinV) / m_sinWminusV);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::YZtoU(const double y, const double z) const
{
    return (z * m_cosU - y * m_sinU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::YZtoV(const double y, const double z) const
{
    return (z * m_cosV - y * m_sinV);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArRotationalTransformationPlugin::YZtoW(const double y, const double z) const
{
    return (z * m_cosW - y * m_sinW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRotationalTransformationPlugin::GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU,
    const double sigmaV, const double sigmaW, double &y, double &z, double &chiSquared) const
{
    const double sigmaU2(sigmaU * sigmaU), sigmaV2(sigmaV * sigmaV), sigmaW2(sigmaW * sigmaW);

    // Obtain expression for chi2, differentiate wrt y and z, set both results to zero and solve simultaneously. Here just paste-in result.
    y = (sigmaW2 * v * m_cosU * m_cosV * m_sinU - sigmaW2 * u * m_cosV * m_cosV * m_sinU + sigmaV2 * w * m_cosU * m_cosW * m_sinU -
            sigmaV2 * u * m_cosW * m_cosW * m_sinU - sigmaW2 * v * m_cosU * m_cosU * m_sinV + sigmaW2 * u * m_cosU * m_cosV * m_sinV +
            sigmaU2 * w * m_cosV * m_cosW * m_sinV - sigmaU2 * v * m_cosW * m_cosW * m_sinV - sigmaV2 * w * m_cosU * m_cosU * m_sinW -
            sigmaU2 * w * m_cosV * m_cosV * m_sinW + sigmaV2 * u * m_cosU * m_cosW * m_sinW + sigmaU2 * v * m_cosV * m_cosW * m_sinW) /
        (sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaV2 * m_cosW * m_cosW * m_sinU * m_sinU - 2. * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV +
            sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaU2 * m_cosW * m_cosW * m_sinV * m_sinV -
            2. * sigmaV2 * m_cosU * m_cosW * m_sinU * m_sinW - 2. * sigmaU2 * m_cosV * m_cosW * m_sinV * m_sinW +
            sigmaV2 * m_cosU * m_cosU * m_sinW * m_sinW + sigmaU2 * m_cosV * m_cosV * m_sinW * m_sinW);

    z = (sigmaW2 * v * m_cosV * m_sinU * m_sinU + sigmaV2 * w * m_cosW * m_sinU * m_sinU - sigmaW2 * v * m_cosU * m_sinU * m_sinV -
            sigmaW2 * u * m_cosV * m_sinU * m_sinV + sigmaW2 * u * m_cosU * m_sinV * m_sinV + sigmaU2 * w * m_cosW * m_sinV * m_sinV -
            sigmaV2 * w * m_cosU * m_sinU * m_sinW - sigmaV2 * u * m_cosW * m_sinU * m_sinW - sigmaU2 * w * m_cosV * m_sinV * m_sinW -
            sigmaU2 * v * m_cosW * m_sinV * m_sinW + sigmaV2 * u * m_cosU * m_sinW * m_sinW + sigmaU2 * v * m_cosV * m_sinW * m_sinW) /
        (sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaV2 * m_cosW * m_cosW * m_sinU * m_sinU - 2. * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV +
            sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaU2 * m_cosW * m_cosW * m_sinV * m_sinV -
            2. * sigmaV2 * m_cosU * m_cosW * m_sinU * m_sinW - 2. * sigmaU2 * m_cosV * m_cosW * m_sinV * m_sinW +
            sigmaV2 * m_cosU * m_cosU * m_sinW * m_sinW + sigmaU2 * m_cosV * m_cosV * m_sinW * m_sinW);

    const double deltaU(u - LArRotationalTransformationPlugin::YZtoU(y, z));
    const double deltaV(v - LArRotationalTransformationPlugin::YZtoV(y, z));
    const double deltaW(w - LArRotationalTransformationPlugin::YZtoW(y, z));
    chiSquared = ((deltaU * deltaU) / sigmaU2) + ((deltaV * deltaV) / sigmaV2) + ((deltaW * deltaW) / sigmaW2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArRotationalTransformationPlugin::GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV,
    const double sigmaW, const double uFit, const double vFit, const double wFit, const double sigmaFit, double &y, double &z, double &chiSquared) const
{
    const double sigmaU2(sigmaU * sigmaU), sigmaV2(sigmaV * sigmaV), sigmaW2(sigmaW * sigmaW), sigmaFit2(sigmaFit * sigmaFit);

    // Obtain expression for chi2, differentiate wrt y and z, set both results to zero and solve simultaneously. Here just paste-in result.
    y = (vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinU + vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU +
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosU * m_cosV * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosU * m_cosV * m_sinU -
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinU - uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU -
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosV * m_cosV * m_sinU - sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosV * m_cosV * m_sinU +
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosW * m_sinU + wFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosW * m_sinU +
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosU * m_cosW * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * w * m_cosU * m_cosW * m_sinU -
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_cosW * m_sinU - uFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_cosW * m_sinU -
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosW * m_cosW * m_sinU - sigmaV2 * sigmaFit2 * sigmaFit2 * u * m_cosW * m_cosW * m_sinU -
            vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinV - vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV -
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosU * m_cosU * m_sinV - sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosU * m_cosU * m_sinV +
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinV + uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinV +
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosU * m_cosV * m_sinV + sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosU * m_cosV * m_sinV +
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosW * m_sinV + wFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosW * m_sinV +
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosV * m_cosW * m_sinV + sigmaU2 * sigmaFit2 * sigmaFit2 * w * m_cosV * m_cosW * m_sinV -
            vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_cosW * m_sinV - vFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_cosW * m_sinV -
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosW * m_cosW * m_sinV - sigmaU2 * sigmaFit2 * sigmaFit2 * v * m_cosW * m_cosW * m_sinV -
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinW - wFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinW -
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosU * m_cosU * m_sinW - sigmaV2 * sigmaFit2 * sigmaFit2 * w * m_cosU * m_cosU * m_sinW -
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinW - wFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinW -
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosV * m_cosV * m_sinW - sigmaU2 * sigmaFit2 * sigmaFit2 * w * m_cosV * m_cosV * m_sinW +
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosW * m_sinW + uFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosU * m_cosW * m_sinW +
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosU * m_cosW * m_sinW + sigmaV2 * sigmaFit2 * sigmaFit2 * u * m_cosU * m_cosW * m_sinW +
            vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosW * m_sinW + vFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosV * m_cosW * m_sinW +
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosV * m_cosW * m_sinW + sigmaU2 * sigmaFit2 * sigmaFit2 * v * m_cosV * m_cosW * m_sinW) /
        (sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
            sigmaV2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
            sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_cosW * m_sinU * m_sinU + sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_cosW * m_sinU * m_sinU +
            sigmaV2 * sigmaW2 * sigmaFit2 * m_cosW * m_cosW * m_sinU * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * m_cosW * m_cosW * m_sinU * m_sinU -
            2. * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV - 2. * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV -
            2. * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV -
            2. * sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV + sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV +
            sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV +
            sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_cosW * m_sinV * m_sinV +
            sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_cosW * m_sinV * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosW * m_cosW * m_sinV * m_sinV +
            sigmaU2 * sigmaFit2 * sigmaFit2 * m_cosW * m_cosW * m_sinV * m_sinV - 2. * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosW * m_sinU * m_sinW -
            2. * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosU * m_cosW * m_sinU * m_sinW - 2. * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosW * m_sinU * m_sinW -
            2. * sigmaV2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosW * m_sinU * m_sinW -
            2. * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosW * m_sinV * m_sinW - 2. * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosV * m_cosW * m_sinV * m_sinW -
            2. * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosW * m_sinV * m_sinW -
            2. * sigmaU2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosW * m_sinV * m_sinW +
            sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinW * m_sinW + sigmaU2 * sigmaV2 * sigmaFit2 * m_cosU * m_cosU * m_sinW * m_sinW +
            sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinW * m_sinW + sigmaV2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosU * m_sinW * m_sinW +
            sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinW * m_sinW + sigmaU2 * sigmaV2 * sigmaFit2 * m_cosV * m_cosV * m_sinW * m_sinW +
            sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinW * m_sinW + sigmaU2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosV * m_sinW * m_sinW);

    z = (vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinU * m_sinU + vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosV * m_sinU * m_sinU +
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosV * m_sinU * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosV * m_sinU * m_sinU +
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_sinU * m_sinU + wFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosW * m_sinU * m_sinU +
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosW * m_sinU * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * w * m_cosW * m_sinU * m_sinU -
            vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinU * m_sinV - vFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_sinU * m_sinV -
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosU * m_sinU * m_sinV - sigmaW2 * sigmaFit2 * sigmaFit2 * v * m_cosU * m_sinU * m_sinV -
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinU * m_sinV - uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_sinU * m_sinV -
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosV * m_sinU * m_sinV - sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosV * m_sinU * m_sinV +
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinV * m_sinV + uFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_sinV * m_sinV +
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosU * m_sinV * m_sinV + sigmaW2 * sigmaFit2 * sigmaFit2 * u * m_cosU * m_sinV * m_sinV +
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_sinV * m_sinV + wFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosW * m_sinV * m_sinV +
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosW * m_sinV * m_sinV + sigmaU2 * sigmaFit2 * sigmaFit2 * w * m_cosW * m_sinV * m_sinV -
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinU * m_sinW - wFit * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_sinU * m_sinW -
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosU * m_sinU * m_sinW - sigmaV2 * sigmaFit2 * sigmaFit2 * w * m_cosU * m_sinU * m_sinW -
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_sinU * m_sinW - uFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_sinU * m_sinW -
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosW * m_sinU * m_sinW - sigmaV2 * sigmaFit2 * sigmaFit2 * u * m_cosW * m_sinU * m_sinW -
            wFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinV * m_sinW - wFit * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_sinV * m_sinW -
            sigmaU2 * sigmaV2 * sigmaFit2 * w * m_cosV * m_sinV * m_sinW - sigmaU2 * sigmaFit2 * sigmaFit2 * w * m_cosV * m_sinV * m_sinW -
            vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_sinV * m_sinW - vFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_sinV * m_sinW -
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosW * m_sinV * m_sinW - sigmaU2 * sigmaFit2 * sigmaFit2 * v * m_cosW * m_sinV * m_sinW +
            uFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_sinW * m_sinW + uFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosU * m_sinW * m_sinW +
            sigmaV2 * sigmaW2 * sigmaFit2 * u * m_cosU * m_sinW * m_sinW + sigmaV2 * sigmaFit2 * sigmaFit2 * u * m_cosU * m_sinW * m_sinW +
            vFit * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_sinW * m_sinW + vFit * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosV * m_sinW * m_sinW +
            sigmaU2 * sigmaW2 * sigmaFit2 * v * m_cosV * m_sinW * m_sinW + sigmaU2 * sigmaFit2 * sigmaFit2 * v * m_cosV * m_sinW * m_sinW) /
        (sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
            sigmaV2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU + sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosV * m_sinU * m_sinU +
            sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_cosW * m_sinU * m_sinU + sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_cosW * m_sinU * m_sinU +
            sigmaV2 * sigmaW2 * sigmaFit2 * m_cosW * m_cosW * m_sinU * m_sinU + sigmaV2 * sigmaFit2 * sigmaFit2 * m_cosW * m_cosW * m_sinU * m_sinU -
            2. * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosV * m_sinU * m_sinV - 2. * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV -
            2. * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV -
            2. * sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosV * m_sinU * m_sinV + sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinV * m_sinV +
            sigmaU2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV +
            sigmaW2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosU * m_sinV * m_sinV + sigmaU2 * sigmaV2 * sigmaW2 * m_cosW * m_cosW * m_sinV * m_sinV +
            sigmaU2 * sigmaV2 * sigmaFit2 * m_cosW * m_cosW * m_sinV * m_sinV + sigmaU2 * sigmaW2 * sigmaFit2 * m_cosW * m_cosW * m_sinV * m_sinV +
            sigmaU2 * sigmaFit2 * sigmaFit2 * m_cosW * m_cosW * m_sinV * m_sinV - 2. * sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosW * m_sinU * m_sinW -
            2. * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosU * m_cosW * m_sinU * m_sinW - 2. * sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosW * m_sinU * m_sinW -
            2. * sigmaV2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosW * m_sinU * m_sinW -
            2. * sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosW * m_sinV * m_sinW - 2. * sigmaU2 * sigmaV2 * sigmaFit2 * m_cosV * m_cosW * m_sinV * m_sinW -
            2. * sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosW * m_sinV * m_sinW -
            2. * sigmaU2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosW * m_sinV * m_sinW +
            sigmaU2 * sigmaV2 * sigmaW2 * m_cosU * m_cosU * m_sinW * m_sinW + sigmaU2 * sigmaV2 * sigmaFit2 * m_cosU * m_cosU * m_sinW * m_sinW +
            sigmaV2 * sigmaW2 * sigmaFit2 * m_cosU * m_cosU * m_sinW * m_sinW + sigmaV2 * sigmaFit2 * sigmaFit2 * m_cosU * m_cosU * m_sinW * m_sinW +
            sigmaU2 * sigmaV2 * sigmaW2 * m_cosV * m_cosV * m_sinW * m_sinW + sigmaU2 * sigmaV2 * sigmaFit2 * m_cosV * m_cosV * m_sinW * m_sinW +
            sigmaU2 * sigmaW2 * sigmaFit2 * m_cosV * m_cosV * m_sinW * m_sinW + sigmaU2 * sigmaFit2 * sigmaFit2 * m_cosV * m_cosV * m_sinW * m_sinW);

    const double outputU(LArRotationalTransformationPlugin::YZtoU(y, z));
    const double outputV(LArRotationalTransformationPlugin::YZtoV(y, z));
    const double outputW(LArRotationalTransformationPlugin::YZtoW(y, z));

    const double deltaU(u - outputU), deltaV(v - outputV), deltaW(w - outputW);
    const double deltaUFit(uFit - outputU), deltaVFit(vFit - outputV), deltaWFit(wFit - outputW);

    chiSquared = ((deltaU * deltaU) / sigmaU2) + ((deltaV * deltaV) / sigmaV2) + ((deltaW * deltaW) / sigmaW2) +
        ((deltaUFit * deltaUFit) / sigmaFit2) + ((deltaVFit * deltaVFit) / sigmaFit2) + ((deltaWFit * deltaWFit) / sigmaFit2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArRotationalTransformationPlugin::Initialize()
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "LArRotationalTransformationPlugin::Initialize - LArTPC description not registered with Pandora as required " << std::endl;
        return STATUS_CODE_NOT_INITIALIZED;
    }

    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    m_thetaU = pFirstLArTPC->GetWireAngleU();
    m_thetaV = pFirstLArTPC->GetWireAngleV();
    m_thetaW = pFirstLArTPC->GetWireAngleW();
    const double sigmaUVW(pFirstLArTPC->GetSigmaUVW());

    m_sinU = std::sin(m_thetaU);
    m_sinV = std::sin(m_thetaV);
    m_sinW = std::sin(m_thetaW);
    m_cosU = std::cos(m_thetaU);
    m_cosV = std::cos(m_thetaV);
    m_cosW = std::cos(m_thetaW);

    m_sinVminusU = std::sin(m_thetaV - m_thetaU);
    m_sinWminusV = std::sin(m_thetaW - m_thetaV);
    m_sinUminusW = std::sin(m_thetaU - m_thetaW);

    if ((std::fabs(m_sinVminusU) < std::numeric_limits<double>::epsilon()) ||
        (std::fabs(m_sinWminusV) < std::numeric_limits<double>::epsilon()) || (std::fabs(m_sinUminusW) < std::numeric_limits<double>::epsilon()))
    {
        std::cout << "LArRotationalTransformationPlugin::Initialize - Equal wire angles; Plugin does not support provided LArTPC configurations. "
                  << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);

        if ((std::fabs(m_thetaU - pLArTPC->GetWireAngleU()) > m_maxAngularDiscrepancyU) ||
            (std::fabs(m_thetaV - pLArTPC->GetWireAngleV()) > m_maxAngularDiscrepancyV) ||
            (std::fabs(m_thetaW - pLArTPC->GetWireAngleW()) > m_maxAngularDiscrepancyW) ||
            (std::fabs(sigmaUVW - pLArTPC->GetSigmaUVW()) > m_maxSigmaDiscrepancy))
        {
            std::cout << "LArRotationalTransformationPlugin::Initialize - Dissimilar drift volumes; Plugin does not support provided LArTPC configurations. "
                      << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LArRotationalTransformationPlugin::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAngularDiscrepancyU", m_maxAngularDiscrepancyU));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAngularDiscrepancyV", m_maxAngularDiscrepancyV));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAngularDiscrepancyW", m_maxAngularDiscrepancyW));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSigmaDiscrepancy", m_maxSigmaDiscrepancy));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
