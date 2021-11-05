/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the gamma start refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArShowerRefinement/GammaStartRefinementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode GammaStartRefinementAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
