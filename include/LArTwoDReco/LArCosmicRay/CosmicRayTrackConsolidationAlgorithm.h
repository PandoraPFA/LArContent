/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray track consolidation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TRACK_CONSOLIDATION_ALGORITHM_H
#define LAR_COSMIC_RAY_TRACK_CONSOLIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  CosmicRayTrackConsolidationAlgorithm class
 */
class CosmicRayTrackConsolidationAlgorithm : public pandora::Algorithm
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

    void SortInputClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &trackClusters,
        pandora::ClusterVector &shortClusters) const;



    void BuildSlidingLinearFits(const pandora::ClusterVector &trackClusters, TwoDSlidingFitResultList &slidingFitResultList) const;



    void GetAssociatedHits(const TwoDSlidingFitResultList &slidingFitResultList, const pandora::ClusterVector &shortClusters,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;


    void GetAssociatedHits(const TwoDSlidingFitResult& slidingFitResult, const pandora::Cluster* pTargetCluster,
        ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const;



    


    pandora::StatusCode RebuildTrackClusters(const ClusterToHitMap &clustersToRebuild) const;

    pandora::StatusCode RebuildShortClusters(const ClusterToHitMap &clustersToRebuild) const;

    pandora::StatusCode BuildNewClusters(const ClusterToHitMap &clustersToRebuild) const;

    pandora::StatusCode BuildNewClusters(const pandora::CaloHitList &inputCaloHitList) const;



    float m_maxClusterLength;

    float m_minTrackLength;

    float m_maxTransverseDisplacement;

    float m_minAssociatedSpan;

    float m_minAssociatedFraction;

    unsigned int m_halfWindowLayers;

    std::string m_clusteringAlgorithmName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackConsolidationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_TRACK_CONSOLIDATION_ALGORITHM_H
