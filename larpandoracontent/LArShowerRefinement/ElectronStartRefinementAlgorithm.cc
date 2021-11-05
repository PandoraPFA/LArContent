/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronStartRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the electron start refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArShowerRefinement/ElectronStartRefinementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ElectronStartRefinementAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronStartRefinementAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
