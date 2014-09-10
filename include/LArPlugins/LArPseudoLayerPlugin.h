/**
 *  @file   LArContent/include/LArPlugins/LArPseudoLayerPlugin.h
 * 
 *  @brief  Header file for the lar pseudo layer plugin class.
 * 
 *  $Log: $
 */
#ifndef LAR_PSEUDO_LAYER_PLUGIN_H
#define LAR_PSEUDO_LAYER_PLUGIN_H 1

#include "Plugins/PseudoLayerPlugin.h"

namespace lar
{

/**
 *  @brief  LArPseudoLayerPlugin class
 */
class LArPseudoLayerPlugin : public pandora::PseudoLayerPlugin
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

#endif // #ifndef LAR_PSEUDO_LAYER_PLUGIN_H
