/**
 *  @file   LArContent/include/LArPseudoLayerCalculator.h
 * 
 *  @brief  Header file for the lar pseudo layer calculator class.
 * 
 *  $Log: $
 */
#ifndef LAR_PSEUDO_LAYER_CALCULATOR_H
#define LAR_PSEUDO_LAYER_CALCULATOR_H 1

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
     *  @brief  Get the z pitch
     * 
     *  @return the z pitch
     */
    virtual float GetZPitch() const = 0;
};

} // namespace lar

#endif // #ifndef LAR_PSEUDO_LAYER_CALCULATOR_H
