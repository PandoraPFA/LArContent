/**
 *  @file   larpandoracontent/LArPlugins/LArPseudoLayerPlugin.cxx
 *
 *  @brief  Implementation of the LAr pseudo layer Plugin class.
 *
 *  $Log: $
 */

#include "Geometry/LArTPC.h"

#include "Helpers/XmlHelper.h"

#include "Managers/GeometryManager.h"

#include "Pandora/Pandora.h"
#include "Pandora/PandoraInputTypes.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"

#include <limits>

namespace lar_content
{

using namespace pandora;

LArPseudoLayerPlugin::LArPseudoLayerPlugin() :
    m_zPitch(std::numeric_limits<float>::max()),
    m_zOffset(0.01f),
    m_zerothLayer(5000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArPseudoLayerPlugin::GetPseudoLayer(const pandora::CartesianVector &positionVector) const
{
    const float zLayer((positionVector.GetZ() + m_zOffset) / m_zPitch + static_cast<float>(m_zerothLayer));

    if (zLayer < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return static_cast<unsigned int>(zLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArPseudoLayerPlugin::Initialize()
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "LArPseudoLayerPlugin::Initialize - LArTPC description not registered with Pandora as required " << std::endl;
        return STATUS_CODE_NOT_INITIALIZED;
    }

    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    m_zPitch = pFirstLArTPC->GetWirePitchW();

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);

        if (std::fabs(m_zPitch - pLArTPC->GetWirePitchW()) > std::numeric_limits<float>::epsilon())
        {
            std::cout << "LArPseudoLayerPlugin::Initialize - Plugin does not support provided LArTPC configurations " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LArPseudoLayerPlugin::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
