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
    m_zPitch = this->GetPandora().GetGeometry()->GetLArTPC().GetWirePitchW();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LArPseudoLayerPlugin::ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
