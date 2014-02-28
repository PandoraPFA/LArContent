/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayTrackConsolidationAlgorithm::Run()
{
    std::cout << "CosmicRayTrackConsolidationAlgorithm::Run()" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
