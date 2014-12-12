/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/CosmicRayBuildingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic-ray building algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/CosmicRayBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CosmicRayBuildingAlgorithm::Run()
{
    std::cout << " --- CosmicRayBuildingAlgorithm::Run() --- " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayBuildingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{   
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
