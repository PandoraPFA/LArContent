/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayShowerMatchingAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMatchingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
