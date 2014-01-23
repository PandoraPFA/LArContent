/**
 *  @file   LArContent/src/LArTwoDSeed/CosmicRaySeedFindingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic-ray seed-finding algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDSeed/CosmicRaySeedFindingAlgorithm.h"

using namespace pandora;

namespace lar
{

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRaySeedFindingAlgorithm::GetSeedClusterList(const ClusterVector &candidateClusters, ClusterList &seedClusterList) const
{
    // TODO: Write the CosmicRaySeedFindingAlgorithm...

    for (ClusterVector::const_iterator iter = candidateClusters.begin(), iterEnd = candidateClusters.end(); iter != iterEnd; ++iter)
	seedClusterList.insert(*iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return SeedFindingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
