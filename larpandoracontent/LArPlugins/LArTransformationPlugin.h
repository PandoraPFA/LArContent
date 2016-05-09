/**
 *  @file   larpandoracontent/LArPlugins/LArTransformationPlugin.h
 * 
 *  @brief  Header file for the transformation plugin interface class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSFORMATION_PLUGIN_H
#define LAR_TRANSFORMATION_PLUGIN_H 1

#include "Pandora/PandoraInputTypes.h"

namespace lar_content
{

/**
 *  @brief  LArTransformationPlugin class
 */
class LArTransformationPlugin
{
public:   

    /** 
     *  @brief  Transform from (U,V) to W position
     *
     *  @param  U the U position
     *  @param  V the V position
     */
     virtual double UVtoW(const double u, const double v) const = 0;

    /** 
     *  @brief  Transform from (V,W) to U position
     *
     *  @param  V the V position
     *  @param  W the W position
     */
     virtual double VWtoU(const double v, const double w) const = 0;

    /** 
     *  @brief  Transform from (W,U) to V position
     *
     *  @param  W the W position
     *  @param  U the U position
     */
     virtual double WUtoV(const double w, const double u) const = 0;

    /** 
     *  @brief  Transform from (U,V) to Y position
     *
     *  @param  U the U position
     *  @param  V the V position
     */
    virtual double UVtoY(const double u, const double v)  const = 0;

    /** 
     *  @brief  Transform from (U,V) to Z position
     *
     *  @param  U the U position
     *  @param  V the V position
     */
     virtual double UVtoZ(const double u, const double v) const = 0;

    /** 
     *  @brief  Transform from (Y,Z) to U position
     * 
     *  @param  Y the Y position
     *  @param  Z the Z position
     */
     virtual double YZtoU(const double y, const double z) const = 0;

    /** 
     *  @brief  Transform from (Y,Z) to V position
     * 
     *  @param  Y the Y position
     *  @param  Z the Z position
     */
     virtual double YZtoV(const double y, const double z) const = 0;

    /** 
     *  @brief  Transform from (pU,pV) to pW direction
     *
     *  @param  pU the pU direction
     *  @param  pV the pV direction
     */
     virtual double PUPVtoPW(const double pu, const double pv) const;

    /** 
     *  @brief  Transform from (pV,pW) to pU direction
     *
     *  @param  pV the pV direction
     *  @param  pW the pW direction
     */
     virtual double PVPWtoPU(const double pv, const double pw) const;

    /** 
     *  @brief  Transform from (pW,pU) to pV direction
     *
     *  @param  pW the pW direction
     *  @param  pU the pU direction
     */
     virtual double PWPUtoPV(const double pw, const double pu) const;

    /** 
     *  @brief  Transform from (pU,pV) to pY position
     *
     *  @param  pU the pU position
     *  @param  pV the pV position
     */
    virtual double PUPVtoPY(const double pu, const double pv)  const;

    /** 
     *  @brief  Transform from (pU,pV) to pZ position
     *
     *  @param  pU the pU position
     *  @param  pV the pV position
     */
     virtual double PUPVtoPZ(const double pu, const double pv) const;

    /** 
     *  @brief  Transform from (pY,pZ) to pU direction
     * 
     *  @param  pU the U component
     *  @param  pV the V component
     */
    virtual double PYPZtoPU(const double py, const double pz) const;

    /** 
     *  @brief  Transform from (pY,pZ) to pV direction
     * 
     *  @param  pU the U component
     *  @param  pV the V component
     */
    virtual double PYPZtoPV(const double py, const double pz) const;

    /** 
     *  @brief  Get resolution, in cm, for calculation of chi2
     * 
     *  @return resolution, in cm, for calculation of chi2
     */
    virtual double GetSigmaUVW() const = 0;

    /** 
     *  @brief  Get the y, z position that yields the minimum chi squared value with respect to specified u, v and w coordinates
     * 
     *  @param  u the u coordinate
     *  @param  v the v coordinate
     *  @param  w the w coordinate
     *  @param  sigmaU the uncertainty in the u coordinate
     *  @param  sigmaV the uncertainty in the v coordinate
     *  @param  sigmaW the uncertainty in the w coordinate
     *  @param  y to receive the y coordinate
     *  @param  z to receive the z coordinate
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV, const double sigmaW,
        double &y, double &z, double &chiSquared) const = 0;

    typedef std::pair<double, pandora::HitType> PositionAndType;

    /** 
     *  @brief  Get the y, z position that corresponds to a projection of two fit positions onto the specific wire associated with a hit
     * 
     *  @param  hitPositionAndType the hit position and hit type
     *  @param  fitPositionAndType1 the first fit position and hit type
     *  @param  fitPositionAndType2 the second fit position and hit type
     *  @param  sigmaHit the uncertainty in the hit coordinate
     *  @param  sigmaFit the uncertainty in the fit coordinates
     *  @param  y to receive the y coordinate
     *  @param  z to receive the z coordinate
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetProjectedYZ(const PositionAndType &hitPositionAndType, const PositionAndType &fitPositionAndType1,
        const PositionAndType &fitPositionAndType2, const double sigmaHit, const double sigmaFit, double &y, double &z, double &chiSquared) const = 0;
    
};

} // namespace lar_content

#endif // #ifndef LAR_TRANSFORMATION_PLUGIN_H
