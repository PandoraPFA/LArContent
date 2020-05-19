/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleSliceSelectionTool.cc
 *
 *  @brief  Implementation of the beam slice selection tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/BeamParticleSliceSelectionTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

BeamParticleSliceSelectionTool::BeamParticleSliceSelectionTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BeamParticleSliceSelectionTool::IsTarget(const MCParticle *const mcParticle) const
{
    return LArMCParticleHelper::IsBeamParticle(mcParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamParticleSliceSelectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return GenericSliceSelectionTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
