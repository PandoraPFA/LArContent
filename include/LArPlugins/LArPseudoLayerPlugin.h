/**
 *  @file   LArContent/LArPlugins/LarPandoraPseudoLayerPlugin.h
 *
 *  @brief  Header file for the lar pseudo layer plugin class.
 *
 *  $Log: $
 */
#ifndef LAR_PSEUDO_LAYER_PLUGIN_H
#define LAR_PSEUDO_LAYER_PLUGIN_H 1

#include "Plugins/PseudoLayerPlugin.h"

namespace lar_content
{

/**
 *  @brief  LarPandoraPseudoLayerPlugin class
 */
class LArPseudoLayerPlugin : public pandora::PseudoLayerPlugin
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  zPitch the wire pitch
     */
    LArPseudoLayerPlugin(const float zPitch);

    /**
     *  @brief  Get pseudolayer corresponding to a specified z coordinate
     *
     *  @param  zCoordinate the z coordinate
     *
     *  @return the pseudolayer
     */
    unsigned int GetPseudoLayer(const float zCoordinate) const;

    /**
     *  @brief  Get the wire pitch
     *
     *  @return the wire pitch
     */
    float GetZPitch() const;

private:
    unsigned int GetPseudoLayer(const pandora::CartesianVector &positionVector) const;
    unsigned int GetPseudoLayerAtIp() const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const float     m_zPitch;       ///< The z pitch
    float           m_zOffset;      ///< The z offset
    unsigned int    m_zerothLayer;  ///< The zeroth layer
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArPseudoLayerPlugin::GetPseudoLayerAtIp() const
{
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPseudoLayerPlugin::GetZPitch() const
{
    return m_zPitch;
}

} // namespace lar

#endif // #ifndef LAR_PSEUDO_LAYER_PLUGIN_H
