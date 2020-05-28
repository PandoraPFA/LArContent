/**
 *  @file   larpandoracontent/LArCheating/CheatingBeamParticleSliceSelectionTool.cc
 *
 *  @brief  Implementation of the beam slice selection tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingBeamParticleSliceSelectionTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

bool CheatingBeamParticleSliceSelectionTool::IsTarget(const MCParticle *const mcParticle) const
{
    return LArMCParticleHelper::IsBeamParticle(mcParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingBeamParticleSliceSelectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return CheatingSliceSelectionTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
