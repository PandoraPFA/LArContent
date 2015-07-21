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
     *  @brief  Constructor
     *
     *  @param  uPitch the wire pitch
     *  @param  vPitch the wire pitch
     *  @param  wPitch the wire pitch
     */
    LArPseudoLayerPlugin(const float uPitch, const float vPitch, const float wPitch);

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
     */
    float GetZPitch() const;    

    /**
     *  @brief  Get the wire pitch
     */
    float GetUPitch() const;    

    /**
     *  @brief  Get the wire pitch
     */
    float GetVPitch() const;    

    /**
     *  @brief  Get the wire pitch
     */
    float GetWPitch() const;

private:
    unsigned int GetPseudoLayer(const pandora::CartesianVector &positionVector) const;
    unsigned int GetPseudoLayerAtIp() const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const float     m_uPitch;       ///< The u pitch
    const float     m_vPitch;       ///< The v pitch
    const float     m_wPitch;       ///< The w pitch

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

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPseudoLayerPlugin::GetUPitch() const
{
    return m_uPitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPseudoLayerPlugin::GetVPitch() const
{
    return m_vPitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPseudoLayerPlugin::GetWPitch() const
{
    return m_wPitch;
}

} // namespace lar

#endif // #ifndef LAR_PSEUDO_LAYER_PLUGIN_H
