/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h
 *
 *  @brief  Header file for the track consolidation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_CONSOLIDATION_ALGORITHM_H
#define LAR_TRACK_CONSOLIDATION_ALGORITHM_H 1

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar_content
{

/**
 *  @brief  TrackConsolidationAlgorithm class
 */
class TrackConsolidationAlgorithm : public TwoDSlidingFitConsolidationAlgorithm
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

    /**
     *  @brief  Default constructor
     */
    TrackConsolidationAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Get the list of hits to be added to track clusters and removed from shower clusters
     *
     *  @param slidingFitResultList  the list of sliding linear fits to track clusters
     *  @param showerClusters  the vector of shower clusters
     *  @param caloHitsToAdd   the output map of hits to be added to clusters
     *  @param caloHitsToRemove   the output map of hits to be removed from clusters
     */
    void GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultList, const pandora::ClusterVector &showerClusters,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;

    /**
     *  @brief Get the list of hits to be added to a track cluster and removed from a shower cluster
     *
     *  @param slidingFitResult  sliding linear fit to track cluster
     *  @param pTargetCluster  shower cluster
     *  @param caloHitsToAdd  the output map of hits to be added to clusters
     *  @param caloHitsToRemove  the output map of hits to be removed from clusters
     */
    void GetReclusteredHits(const TwoDSlidingFitResult& slidingFitResult, const pandora::Cluster* pTargetCluster,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;

    float        m_maxTransverseDisplacement;  ///<
    float        m_minAssociatedSpan;          ///<
    float        m_minAssociatedFraction;      ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TrackConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new TrackConsolidationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_CONSOLIDATION_ALGORITHM_H
