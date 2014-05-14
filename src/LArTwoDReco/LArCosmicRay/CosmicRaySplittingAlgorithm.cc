/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRaySplittingAlgorithm::Run()
{
   std::cout << " --- CosmicRaySplittingAlgorithm::Run() --- " << std::endl;

   return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
