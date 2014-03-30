/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/TrackCleaningAlgorithm.h
 *
 *  @brief  Header file for the track cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_CLEANING_ALGORITHM_H
#define LAR_TRACK_CLEANING_ALGORITHM_H 1

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  TrackCleaningAlgorithm class
 */
class TrackCleaningAlgorithm : public TwoDSlidingFitConsolidationAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Get the list of hits to be added or removed from clusters
     *
     *  @param slidingFitResultList  the list of sliding linear fits to track clusters
     *  @param showerClusters  the vector of shower clusters
     *  @param caloHitsToAdd   the output map of hits to be added to clusters
     *  @param caloHitsToRemove   the output map of hits to be removed from clusters
     */
    void GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultList, const pandora::ClusterVector &showerClusters,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TrackCleaningAlgorithm::Factory::CreateAlgorithm() const
{
    return new TrackCleaningAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRACK_CLEANING_ALGORITHM_H
