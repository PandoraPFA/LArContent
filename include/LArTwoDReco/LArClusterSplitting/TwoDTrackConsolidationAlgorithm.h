/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/TwoDTrackConsolidationAlgorithm.h
 *
 *  @brief  Header file for the 2D track consolidation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_TRACK_CONSOLIDATION_ALGORITHM_H
#define LAR_TWO_D_TRACK_CONSOLIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  TwoDTrackConsolidationAlgorithm class
 */
class TwoDTrackConsolidationAlgorithm : public pandora::Algorithm
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const pandora::Cluster*, pandora::CaloHitList> ClusterToHitMap;
    typedef std::map<const pandora::CaloHit*, pandora::CaloHitList> HitToHitMap;

    /**
     *  @brief Sort input cluster list into track clusters and shower clusters
     *
     *  @param pClusterList  the input cluster list
     *  @param trackClusters  the output vector of track clusters
     *  @param showerClusters  the output vector of shower clusters
     */
    void SortInputClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &trackClusters,
        pandora::ClusterVector &showerClusters) const;

    /**
     *  @brief Apply sliding linear fits to track clusters
     *
     *  @param trackClusters  the input vector of track clusters
     *  @param slidingFitResultList  the output list of sliding linear fits
     */
    void BuildSlidingLinearFits(const pandora::ClusterVector &trackClusters, TwoDSlidingFitResultList &slidingFitResultList) const;

    /**
     *  @brief Get the list of hits to be removed from shower clusters and added to track clusters
     *
     *  @param slidingFitResultList  the list of sliding linear fits to track clusters
     *  @param showerClusters  the vector of shower clusters
     *  @param caloHitsToAdd   the output map of hits to be added to the track clusters
     *  @param caloHitsToRemove   the output map of hits to be removed from the shower clusters
     */
    void GetAssociatedHits(const TwoDSlidingFitResultList &slidingFitResultList, const pandora::ClusterVector &showerClusters,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;

    /**
     *  @brief Get the list of hits to be removed from a shower cluster and added to a track cluster
     *
     *  @param slidingFitResult  sliding linear fit to track cluster
     *  @param pTargetCluster  shower cluster
     *  @param caloHitsToAdd  the output map of hits to be added to the track clusters
     *  @param caloHitsToRemove  the output map of hits to be removed from the shower clusters
     */
    void GetAssociatedHits(const TwoDSlidingFitResult& slidingFitResult, const pandora::Cluster* pTargetCluster,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;

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

    /**
     *  @brief Remove hits from shower clusters
     *
     *  @param clustersToRebuild  the list of hits to be removed from each cluster
     *  @param modifiedShowers  the output list of modified showers
     */
    pandora::StatusCode RemoveHitsFromShowers(const ClusterToHitMap &clustersToRebuild, pandora::ClusterList &modifiedShowers) const;

    /**
     *  @brief Add hits to track clusters
     *
     *  @param clustersToRebuild  the list of hits to be added to each cluster
     *  @param modifiedTracks  the output list of modified tracks
     */
    pandora::StatusCode AddHitsToTracks(const ClusterToHitMap &clustersToRebuild, pandora::ClusterList &modifiedTracks) const;

    /**
     *  @brief Re-build clusters
     *
     *  @param tracksToRebuild  the list of tracks to rebuild
     *  @param showersToRebuild  the list of showers to rebuild
     */
    pandora::StatusCode RebuildClusters(const pandora::ClusterList &tracksToRebuild, const pandora::ClusterList &showersToRebuild) const;

    /**
     *  @brief Re-build shower clusters
     *
     *  @param clustersToRebuild  the list of hits to be kept for each cluster
     */
    pandora::StatusCode BuildNewClusters(const ClusterToHitMap &clustersToRebuild) const;

    /**
     *  @brief Re-build clusters
     *
     *  @param inputCaloHitList  the input list of hits
     */
    pandora::StatusCode BuildNewClusters(const pandora::CaloHitList &inputCaloHitList) const;


    float        m_maxClusterLength;           ///<
    float        m_minTrackLength;             ///<
    float        m_maxTransverseDisplacement;  ///<
    float        m_minAssociatedSpan;          ///<
    float        m_minAssociatedFraction;      ///<
    float        m_reclusteringWindow;         ///<
    unsigned int m_halfWindowLayers;           ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TwoDTrackConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new TwoDTrackConsolidationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TWO_D_TRACK_CONSOLIDATION_ALGORITHM_H
