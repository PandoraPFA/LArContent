/**
 *  @file   LArContent/include/LArBFieldCalculator.h
 * 
 *  @brief  Header file for the lar bfield calculator class.
 * 
 *  $Log: $
 */

#ifndef LAR_BFIELD_CALCULATOR_H
#define LAR_BFIELD_CALCULATOR_H 1

#include "Utilities/BFieldCalculator.h"

namespace lar
{

/**
 *  @brief  LArBFieldCalculator class
 */
class LArBFieldCalculator : public pandora::BFieldCalculator
{
private:
    void Initialize(const pandora::GeometryHelper *const pGeometryHelper);
    float GetBField(const pandora::CartesianVector &positionVector) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArBFieldCalculator::Initialize(const pandora::GeometryHelper *const pGeometryHelper)
{
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArBFieldCalculator::GetBField(const pandora::CartesianVector &/*positionVector*/) const
{
    return 0.f;
};

} // namespace lar

#endif // #ifndef LAR_BFIELD_CALCULATOR_H
