/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the shower characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArShowerRefinement/ShowerCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ShowerCharacterisationAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCharacterisationAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
