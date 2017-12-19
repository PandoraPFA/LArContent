/**
 *  @file   larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h
 *
 *  @brief  Header file for the rotational transformation plugin class.
 *
 *  $Log: $
 */
#ifndef LAR_ROTATIONAL_TRANSFORMATION_PLUGIN_H
#define LAR_ROTATIONAL_TRANSFORMATION_PLUGIN_H 1

#include "Plugins/LArTransformationPlugin.h"

namespace lar_content
{

/**
 *  @brief  LArRotationalTransformationPlugin class
 */
class LArRotationalTransformationPlugin : public pandora::LArTransformationPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    LArRotationalTransformationPlugin();

    virtual double UVtoW(const double u, const double v) const;
    virtual double VWtoU(const double v, const double w) const;
    virtual double WUtoV(const double w, const double u) const;
    virtual double UVtoY(const double u, const double v) const;
    virtual double UVtoZ(const double u, const double v) const;
    virtual double YZtoU(const double y, const double z) const;
    virtual double YZtoV(const double y, const double z) const;
    virtual double GetSigmaUVW() const;
    virtual void GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV, const double sigmaW,
        double &y, double &z, double &chiSquared) const;
    virtual void GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV, const double sigmaW,
        const double uFit, const double vFit, const double wFit, const double sigmaFit, double &y, double &z, double &chiSquared) const;
    virtual void GetProjectedYZ(const PositionAndType &hitPositionAndType, const PositionAndType &fitPositionAndType1,
        const PositionAndType &fitPositionAndType2, const double sigmaHit, const double sigmaFit, double &y, double &z, double &chiSquared) const;

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    double    m_thetaU;          ///< inclination of U wires (radians)
    double    m_thetaV;          ///< inclination of V wires (radians)
    double    m_sigmaUVW;        ///< resolution (cm), for calculation of chi2
    double    m_sinUplusV;       ///< sin(thetaU+thetaV)
    double    m_sinU;            ///< sin(thetaU)
    double    m_sinV;            ///< sin(thetaV)
    double    m_cosU;            ///< cos(thetaU)
    double    m_cosV;            ///< cos(thetaV)

    double    m_maxAngularDiscrepancy;      ///< Maximum allowed difference between like wire (u and v) angles between LArTPCs
    double    m_maxSigmaDiscrepancy;        ///< Maximum allowed difference between like wire sigma values between LArTPCs
};

} // namespace lar_content

#endif // #ifndef LAR_ROTATIONAL_TRANSFORMATION_PLUGIN_H
