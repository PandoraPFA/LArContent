/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TrackCleaningAlgorithm.cc
 *
 *  @brief  Implementation of the track cleaning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/TrackCleaningAlgorithm.h"

using namespace pandora;

namespace lar
{

void TrackCleaningAlgorithm::GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultListI,
    const ClusterVector &showerClustersJ, ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
   
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackCleaningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
   

    return TwoDSlidingFitConsolidationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
