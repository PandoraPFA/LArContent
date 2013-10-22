/**
 *  @file   LArContent/include/LArHelpers/LArTransformationCalculator.h
 * 
 *  @brief  Header file for the transformation calculator interface class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSFORMATION_CALCULATOR_H
#define LAR_TRANSFORMATION_CALCULATOR_H 1

namespace lar
{

/**
 *  @brief  LArTransformationCalculator class
 */
class LArTransformationCalculator
{
public:
    /** 
     *  @brief  Transform from (U,V) to W position
     *
     *  @param  U the U position
     *  @param  V the V position
     */
     virtual float UVtoW(const float u, const float v) const = 0;

    /** 
     *  @brief  Transform from (V,W) to U position
     *
     *  @param  V the V position
     *  @param  W the W position
     */
     virtual float VWtoU(const float v, const float w) const = 0;

    /** 
     *  @brief  Transform from (W,U) to V position
     *
     *  @param  W the W position
     *  @param  U the U position
     */
     virtual float WUtoV(const float w, const float u) const = 0;

    /** 
     *  @brief  Transform from (U,V) to Y position
     *
     *  @param  U the U position
     *  @param  V the V position
     */
    virtual float UVtoY(const float u, const float v)  const = 0;

    /** 
     *  @brief  Transform from (U,V) to Z position
     *
     *  @param  U the U position
     *  @param  V the V position
     */
     virtual float UVtoZ(const float u, const float v) const = 0;

    /** 
     *  @brief  Transform from (Y,Z) to U position
     * 
     *  @param  Y the Y position
     *  @param  Z the Z position
     */
     virtual float YZtoU(const float y, const float z) const = 0;

    /** 
     *  @brief  Transform from (Y,Z) to V position
     * 
     *  @param  Y the Y position
     *  @param  Z the Z position
     */
     virtual float YZtoV(const float y, const float z) const = 0;

    /** 
     *  @brief  Transform from (pU,pV) to pW direction
     *
     *  @param  pU the pU direction
     *  @param  pV the pV direction
     */
     virtual float PUPVtoPW(const float pu, const float pv) const = 0;

    /** 
     *  @brief  Transform from (pV,pW) to pU direction
     *
     *  @param  pV the pV direction
     *  @param  pW the pW direction
     */
     virtual float PVPWtoPU(const float pv, const float pw) const = 0;

    /** 
     *  @brief  Transform from (pW,pU) to pV direction
     *
     *  @param  pW the pW direction
     *  @param  pU the pU direction
     */
     virtual float PWPUtoPV(const float pw, const float pu) const = 0;

    /** 
     *  @brief  Transform from (pY,pZ) to pU direction
     * 
     *  @param  pU the U component
     *  @param  pV the V component
     */
    virtual float PYPZtoPU(const float py, const float pz) const = 0;

    /** 
     *  @brief  Transform from (pY,pZ) to pV direction
     * 
     *  @param  pU the U component
     *  @param  pV the V component
     */
    virtual float PYPZtoPV(const float py, const float pz) const = 0;

    /** 
     *  @brief  Get resolution, in cm, for calculation of chi2
     * 
     *  @return resolution, in cm, for calculation of chi2
     */
    virtual float GetSigmaUVW() const = 0;
};

} // namespace lar

#endif // #ifndef LAR_TRANSFORMATION_CALCULATOR_H
