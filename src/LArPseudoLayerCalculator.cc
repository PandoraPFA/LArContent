/**
 *  @file   LArPseudoLayerCalculator.cxx
 * 
 *  @brief  Implementation of the lar pseudo layer calculator class.
 * 
 *  $Log: $
 */

#include "LArPseudoLayerCalculator.h"

namespace lar
{

const float LArPseudoLayerCalculator::Z_PITCH = 0.3f;
const float LArPseudoLayerCalculator::Z_OFFSET = 0.01f;

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPseudoLayerCalculator::GetZCoordinate(const pandora::PseudoLayer pseudoLayer)
{
    if (0 == pseudoLayer)
        return 0.f;

    const float zCoordinate((static_cast<float>(pseudoLayer - 1) * Z_PITCH));

    return zCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::PseudoLayer LArPseudoLayerCalculator::GetPseudoLayer(const float zCoordinate)
{
    if (zCoordinate < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const pandora::PseudoLayer pseudoLayer(static_cast<unsigned int>((zCoordinate + Z_OFFSET) / Z_PITCH) + 1);

    return pseudoLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPseudoLayerCalculator::Initialize(const pandora::GeometryHelper *const pGeometryHelper)
{
    // No initialization required
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::PseudoLayer LArPseudoLayerCalculator::GetPseudoLayer(const pandora::CartesianVector &positionVector) const
{
    return LArPseudoLayerCalculator::GetPseudoLayer(positionVector.GetZ());
}

} // namespace lar
