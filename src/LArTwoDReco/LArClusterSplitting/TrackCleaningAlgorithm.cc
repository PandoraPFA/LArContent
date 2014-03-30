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

void TrackCleaningAlgorithm::GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultList,
    const ClusterVector &showerClusters, ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const
{
    // TODO: WRITE THIS METHOD!
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackCleaningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return TwoDSlidingFitConsolidationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
