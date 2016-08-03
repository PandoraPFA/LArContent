/**
 *  @file   larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h
 *
 *  @brief  Header file for the rotational transformation plugin class.
 *
 *  $Log: $
 */
#ifndef LAR_ROTATIONAL_TRANSFORMATION_PLUGIN_H
#define LAR_ROTATIONAL_TRANSFORMATION_PLUGIN_H 1

#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

namespace lar_content
{

/**
 *  @brief  LArRotationalTransformationPlugin class
 */
class LArRotationalTransformationPlugin : public lar_content::LArTransformationPlugin
{
public:

    /**
     *  @brief  Constructor
     *
     *  @param  thetaU  the angle from the W axis to the U axis (radians)
     *  @param  thetaV  the angle from the W axis to the V axis (radians)
     *  @param  sigmaUVW  nominal spatial resolution for U, V and W (cm)
     */
    LArRotationalTransformationPlugin(const double thetaU, const double thetaV, const double sigmaUVW);

    /**
     *  @brief  Destructor
     */
    virtual ~LArRotationalTransformationPlugin();

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

    virtual void GetProjectedYZ(const PositionAndType &hitPositionAndType, const PositionAndType &fitPositionAndType1,
        const PositionAndType &fitPositionAndType2, const double sigmaHit, const double sigmaFit, double &y, double &z, double &chiSquared) const;

private:

    const double    m_thetaU;          ///< inclination of U wires (radians)
    const double    m_thetaV;          ///< inclination of V wires (radians)
    const double    m_sigmaUVW;        ///< resolution (cm), for calculation of chi2
    const double    m_sinUplusV;       ///< sin(thetaU+thetaV)
    const double    m_sinU;            ///< sin(thetaU)
    const double    m_sinV;            ///< sin(thetaV)
    const double    m_cosU;            ///< cos(thetaU)
    const double    m_cosV;            ///< cos(thetaV)
};

} // namespace lar_content

#endif // #ifndef LAR_ROTATIONAL_TRANSFORMATION_PLUGIN_H
