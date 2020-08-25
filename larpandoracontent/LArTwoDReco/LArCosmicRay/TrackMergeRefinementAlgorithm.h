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
class TrackMergeRefinementAlgorithm :  public TrackRefinementBaseAlgorithm
{
public:
   /**
     *  @brief  Default constructor
     */        
    TrackMergeRefinementAlgorithm();
    
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find the best cluster association
     *
     *  @param  clusterVector the vector of clusters to consider
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     *  @param  clusterAssociation the cluster pair association
     *
     *  @return  bool whether a cluster pair association was found
     */    
    bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterPairAssociation &clusterAssociation) const;

    /**
     *  @brief  Whether two clusters are assoicated to one another
     *
     *  @param  upstreamPoint the merge point of the upstream cluster
     *  @param  upstreamDirection the local direction of the upstream cluster at the merge point
     *  @param  downstreamPoint the merge point of the downstream cluster
     *  @param  downstreamDirection the local direction of the downstrean cluster at the merge point
     *
     *  @return  bool whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::CartesianVector &upstreamPoint, const pandora::CartesianVector &upstreamDirection,
        const pandora::CartesianVector &downstreamPoint, const pandora::CartesianVector &downstreamDirection) const;

    /**
     *  @brief  Collect the hits that lie close to the line connecting the associated clusters
     *
     *  @param  clusterAssociation the cluster pair association
     *  @param  pClusterList the list of all clusters
     *  @param  createdMainTrackClusters the list of previously created main track clusters whose hits cannot be collected by other clusters
     *  @param  clusterToCaloHitListMap the output map [parent cluster -> list of hits which belong to the main track]
     */    
    void GetExtrapolatedCaloHits(ClusterPairAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList,
        const pandora::ClusterList &createdMainTrackClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const;
    
    /**
     *  @brief  Perform topological checks on the collected hits to ensure no gaps are present
     *
     *  @param  clusterAssociation the cluster pair association
     *  @param  clusterToCaloHitListMap the input map [parent cluster -> list of hits which belong to the main track]
     *
     *  @return  bool wether the checks pass
     */        
    bool AreExtrapolatedHitsGood(const ClusterPairAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap) const;
  
    /**
     *  @brief  Check the separation of the extremal extrapolated hits with the cluster merge points or, in the case of no hits, the cluster merge point separation
     *
     *  @param  clusterAssociation the cluster pair association
     *  @param  clusterHitVector the extrapolated hit vector (ordered closest hit to the upstream merge point -> furthest hit)
     *
     *  @return  bool wether the checks pass
     */          
    bool AreExtrapolatedHitsNearMergePoints(const ClusterPairAssociation &clusterAssociation, const pandora::CaloHitVector &extrapolatedHitVector) const;

    /**
     *  @brief  Remove the cluster association from the cluster vector so that the same cluster pair is not considered again
     *
     *  @param  clusterAssociation the cluster pair association
     *  @param  clusterVector the vector of clusters considered in future iterations of the algorithm
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     */           
    void ConsiderClusterAssociation(const ClusterPairAssociation &clusterAssociation, pandora::ClusterVector &clusterVector,
        SlidingFitResultMapPair &slidingFitResultMapPair) const;

    /**
     *  @brief  Refine the cluster endpoints and merge together the associated clusters alongside any extrapolated hits
     *
     *  @param  clusterAssociation the cluster pair association
     *  @param  clusterToCaloHitListMap the map [parent cluster -> list of hits which belong to the main track]
     *  @param  pClusterList the list of all clusters
     *  @param  clusterVector the vector of clusters considered in future iterations of the algorithm
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     *
     *  @return  Cluster the address of the created main track cluster
     */    
    const pandora::Cluster *CreateMainTrack(const ClusterPairAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
        const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    unsigned int m_maxLoopIterations;         ///< The maximum number of main loop iterations
    float m_minClusterLengthSum;              ///< The threshold cluster and associated cluster length sum
    float m_minSeparationDistance;            ///< The threshold separation distance between associated clusters
    float m_minDirectionDeviationCosAngle;    ///< The threshold cos opening angle of the associated cluster directions
    float m_maxPredictedMergePointOffset;     ///< The threshold separation distance between the predicted and true cluster merge points
    float m_distanceFromLine;                 ///< The threshold hit distance of an extrapolated hit from the cluster connecting line
    float m_boundaryTolerance;                ///< The maximum allowed distance of an extremal extrapolate hit to a cluster merge point

};
    
} // namespace lar_content

#endif // #ifndef LAR_TRACK_MERGE_REFINEMENT_ALGORITHM_H
