/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoSliceSelectionTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoSliceSelectionTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

bool CheatingNeutrinoSliceSelectionTool::IsTarget(const MCParticle *const mcParticle) const
{
    return LArMCParticleHelper::IsNeutrino(mcParticle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoSliceSelectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return CheatingSliceSelectionTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
