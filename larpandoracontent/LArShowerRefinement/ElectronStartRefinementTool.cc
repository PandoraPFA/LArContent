/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.cc
 *
 *  @brief  Implementation of the electron start refinement tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

ElectronStartRefinementTool::ElectronStartRefinementTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
    const pandora::CartesianVector &nuVertexPosition)
{
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronStartRefinementTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
