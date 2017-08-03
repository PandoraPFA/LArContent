/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/CosmicRayTaggingTool.cc
 *
 *  @brief  Implementation of the cosmic-ray tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/CosmicRayTaggingTool.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTaggingTool::CosmicRayTaggingTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &/*parentCosmicRayPfos*/, PfoList &/*ambiguousPfos*/) const
{
    // TODO
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTaggingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
