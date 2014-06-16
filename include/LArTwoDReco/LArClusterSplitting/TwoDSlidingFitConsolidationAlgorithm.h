/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.h
 *
 *  @brief  Header file for the 2D sliding fit consolidation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_CONSOLIDATION_ALGORITHM_H
#define LAR_TWO_D_SLIDING_FIT_CONSOLIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  TwoDSlidingFitConsolidationAlgorithm class
 */
class TwoDSlidingFitConsolidationAlgorithm : public pandora::Algorithm
{
protected:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const pandora::Cluster*, pandora::CaloHitList> ClusterToHitMap;
    typedef std::map<const pandora::CaloHit*, pandora::CaloHitList> HitToHitMap;

    /**
     *  @brief Get the list of hits to be added or removed from clusters
     *
     *  @param slidingFitResultList  the list of sliding linear fits to track clusters
     *  @param showerClusters  the vector of shower clusters
     *  @param caloHitsToAdd   the output map of hits to be added to clusters
     *  @param caloHitsToRemove   the output map of hits to be removed from clusters
     */
    virtual void GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultList, const pandora::ClusterVector &showerClusters,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const = 0;

private:

    /**
     *  @brief Sort input cluster list into track-like clusters and shower-like clusters
     *
     *  @param pClusterList  the input cluster list
     *  @param trackClusters  the output vector of track-like clusters
     *  @param showerClusters  the output vector of shower-like clusters
     */
    void SortInputClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &trackClusters,
        pandora::ClusterVector &showerClusters) const;

    /**
     *  @brief Apply sliding linear fits to track clusters
     *
     *  @param trackClusters  the input vector of track-like clusters
     *  @param slidingFitResultList  the output list of sliding linear fits
     */
    void BuildSlidingLinearFits(const pandora::ClusterVector &trackClusters, TwoDSlidingFitResultList &slidingFitResultList) const;

    /**
     *  @brief Remove hits from clusters
     *
     *  @param clustersToRebuild  the list of hits to be removed from clusters
     *  @param unavailableClusters  the list of deleted clusters
     */
    pandora::StatusCode RemoveHitsFromClusters(const ClusterToHitMap &clustersToRebuild, pandora::ClusterList &unavailableClusters) const;

    /**
     *  @brief Add hits to clusters
     *
     *  @param clustersToRebuild  the list of hits to be added to clusters
     *  @param unavailableClusters the list of modified clusters
     */
    pandora::StatusCode AddHitsToClusters(const ClusterToHitMap &clustersToRebuild, pandora::ClusterList &unavailableClusters) const;

    /**
     *  @brief Re-build clusters
     *
     *  @param clustersAtStart  the initial mapping of clusters to hits
     *  @param unavailableClusters  the list of unavailable clusters
     */
    pandora::StatusCode RebuildClusters(const ClusterToHitMap &clustersAtStart, const pandora::ClusterList &unavailableClusters) const;


    std::string  m_reclusteringAlgorithmName;  ///< Name of daughter algorithm to use for cluster re-building
    float        m_maxClusterLength;           ///< Maximum length of shower clusters to use in re-building
    float        m_minTrackLength;             ///< Minimum length of track clusters to consolidate
    unsigned int m_halfWindowLayers;           ///< Size of layer window for sliding fit results
};

} // namespace lar

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_CONSOLIDATION_ALGORITHM_H
