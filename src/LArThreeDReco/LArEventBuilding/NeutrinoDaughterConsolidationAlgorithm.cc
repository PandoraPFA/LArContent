/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/NeutrinoDaughterConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino daughter consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/NeutrinoDaughterConsolidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoDaughterConsolidationAlgorithm::NeutrinoDaughterConsolidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoDaughterConsolidationAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoDaughterConsolidationAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
