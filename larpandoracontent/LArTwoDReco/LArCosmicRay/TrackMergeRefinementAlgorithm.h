/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackMergeRefinementAlgorithm.h
 *
 *  @brief  Header file for the track merge refinement class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_MERGE_REFINEMENT_ALGORITHM_H
#define LAR_TRACK_MERGE_REFINEMENT_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackRefinementBaseAlgorithm.h"

namespace lar_content
{
/**
 *  @brief TrackMergeRefinementAlgorithm class
 */
class TrackMergeRefinementAlgorithm :  public TrackRefinementBaseAlgorithm<ClusterPairAssociation>
{
public:
    
    TrackMergeRefinementAlgorithm();
    
protected:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterPairAssociation &clusterAssociation) const;
    
    bool AreClustersAssociated(const pandora::CartesianVector &upstreamPoint, const pandora::CartesianVector &upstreamDirection, const pandora::CartesianVector &downstreamPoint,
        const pandora::CartesianVector &downstreamDirection) const;
    
    void GetExtrapolatedCaloHits(ClusterPairAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, const pandora::ClusterList &createdMainTrackClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    //const pandora::Cluster *CreateMainTrack(ClusterPairAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    void ConsiderClusterAssociation(const ClusterPairAssociation &clusterAssociation, pandora::ClusterVector &clusterVector, pandora::ClusterList &consideredClusters,
        SlidingFitResultMapPair &slidingFitResultMapPair) const;

    //bool AreExtrapolatedHitsGood(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const bool isHigherXBoundary) const;

    float m_minClusterLengthSum;
    float m_minSeparationDistance;
    float m_minDirectionDeviationCosAngle;
    float m_maxPredictedMergePointOffset;    
};
    
} // namespace lar_content

#endif // #ifndef LAR_TRACK_MERGE_REFINEMENT_ALGORITHM_H
