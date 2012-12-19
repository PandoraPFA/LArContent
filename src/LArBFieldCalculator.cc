/**
 *  @file   LArBFieldCalculator.cc
 * 
 *  @brief  Implementation of the lar bfield calculator class.
 * 
 *  $Log: $
 */

#include "LArBFieldCalculator.h"

namespace lar
{

void LArBFieldCalculator::Initialize(const pandora::GeometryHelper *const pGeometryHelper)
{
};

//------------------------------------------------------------------------------------------------------------------------------------------

float LArBFieldCalculator::GetBField(const pandora::CartesianVector &positionVector) const
{
    return 0.f;
};

} // namespace lar
