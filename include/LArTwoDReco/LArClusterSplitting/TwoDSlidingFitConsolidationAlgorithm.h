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
     *  @brief Get initial hit configuration
     *
     *  @param inputClusters  the input vector of clusters
     *  @param caloHitsToStart  the output map of clusters to hits
     */
    void GetInitialHits(const pandora::ClusterVector &inputClusters, ClusterToHitMap &caloHitsAtStart) const;

    /**
     *  @brief Remove hits from clusters
     *
     *  @param clustersToRebuild  the list of hits to be removed from clusters
     */
    pandora::StatusCode RemoveHitsFromClusters(const ClusterToHitMap &clustersToRebuild) const;

    /**
     *  @brief Add hits to clusters
     *
     *  @param clustersToRebuild  the list of hits to be added to clusters
     */
    pandora::StatusCode AddHitsToClusters(const ClusterToHitMap &clustersToRebuild) const;

    /**
     *  @brief Re-build clusters
     *
     *  @param clustersAtStart  the initial mapping of clusters to hits
     */
    pandora::StatusCode RebuildClusters(const ClusterToHitMap &clustersAtStart) const;

    /**
     *  @brief Re-build shower clusters
     *
     *  @param clustersToRebuild  the mapping of initial clusters to their available hits
     */
    pandora::StatusCode BuildNewClusters(const ClusterToHitMap &clustersToRebuild) const;

    /**
     *  @brief Re-build clusters
     *
     *  @param inputCaloHitList  the input list of available hits
     */
    pandora::StatusCode BuildNewClusters(const pandora::CaloHitList &inputCaloHitList) const;

    /**
     *  @brief Re-clustering algorithm used to re-build shower clusters
     *
     *  @param pSeedCaloHit  the seed hit
     *  @param pCurrentCaloHit  a candidate hit
     *  @param hitAssociationMap  the map of associations between all hits
     *  @param vetoList  the list of used hits
     *  @param mergeList  the list of hits associated with the seed hit
     */
    void CollectAssociatedHits(pandora::CaloHit *pSeedCaloHit, pandora::CaloHit *pCurrentCaloHit,
        const HitToHitMap &hitAssociationMap, const pandora::CaloHitList &vetoList, pandora::CaloHitList &mergeList) const;

    float        m_maxClusterLength;           ///<
    float        m_minTrackLength;             ///<
    float        m_reclusteringWindow;         ///<
    unsigned int m_halfWindowLayers;           ///<
};

} // namespace lar

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_CONSOLIDATION_ALGORITHM_H
