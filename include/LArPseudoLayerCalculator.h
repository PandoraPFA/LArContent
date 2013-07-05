/**
 *  @file   LArContent/include/LArPseudoLayerCalculator.h
 * 
 *  @brief  Header file for the lar pseudo layer calculator class.
 * 
 *  $Log: $
 */
#ifndef LAR_PSEUDO_LAYER_CALCULATOR_H
#define LAR_PSEUDO_LAYER_CALCULATOR_H 1

#include "Helpers/GeometryHelper.h"
#include "Utilities/PseudoLayerCalculator.h"

namespace lar
{

/**
 *  @brief  LArPseudoLayerCalculator class
 */
class LArPseudoLayerCalculator : public pandora::PseudoLayerCalculator
{
public:
    /**
     *  @brief  Get z coordinate corresponding to a specified pseudolayer
     * 
     *  @param  pseudolayer the pseudolayer
     * 
     *  @return the z coordinate
     */
    static float GetZCoordinate(const pandora::PseudoLayer pseudoLayer);

    /**
     *  @brief  Get pseudolayer corresponding to a specified z coordinate
     * 
     *  @param  zCoordinate the z coordinate
     * 
     *  @return the pseudolayer
     */
    static pandora::PseudoLayer GetPseudoLayer(const float zCoordinate);

private:
    void Initialize(const pandora::GeometryHelper *const pGeometryHelper);
    pandora::PseudoLayer GetPseudoLayer(const pandora::CartesianVector &positionVector) const;
    pandora::PseudoLayer GetPseudoLayerAtIp() const;

    static const float  Z_PITCH;    ///< The z pitch
    static const float  Z_OFFSET;   ///< The z offset
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::PseudoLayer LArPseudoLayerCalculator::GetPseudoLayerAtIp() const
{
    return 0;
}

} // namespace lar

#endif // #ifndef LAR_PSEUDO_LAYER_CALCULATOR_H
