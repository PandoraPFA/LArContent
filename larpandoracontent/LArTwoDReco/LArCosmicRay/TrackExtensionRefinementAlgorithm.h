/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackExtensionRefinementAlgorithm.h
 *
 *  @brief  Header file for the track extension refinement class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_EXTENSION_REFINEMENT_ALGORITHM_H
#define LAR_TRACK_EXTENSION_REFINEMENT_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackRefinementBaseAlgorithm.h"

namespace lar_content
{
/**
 *  @brief TrackExtensionRefinementAlgorithm class
 */
class TrackExtensionRefinementAlgorithm :  public TrackRefinementBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */    
    TrackExtensionRefinementAlgorithm();
    
protected:
    pandora::StatusCode Run();    
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;

    /**
     *  @brief  Initialise the detector and TPC edge member variables taking into account the gap in the case of a CPA
     */       
    void InitialiseGeometry();

    /**
     *  @brief  Find the best cluster endpoint association
     *
     *  @param  clusterVector the vector of clusters to consider
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     *  @param  pClusterList the list of all clusters
     *  @param  isHigherXboundary whether to look for endpoints that are closer to the higher or lower x tpc boundary
     *  @param  clusterAssociation the output cluster endpoint association
     *
     *  @return bool whether a cluster endpoint association was found
     */    
    virtual bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        const pandora::ClusterList *const pClusterList, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation) = 0;    

    /**
     *  @brief  Collect the extrapolated hits beyond the cluster endpoint that lie close to the projected path of the cluster
     *
     *  @param  clusterAssociation the cluster endpoint association
     *  @param  pClusterList the list of all clusters
     *  @param  createdMainTrackClusters the list of previously created main track clusters whose hits cannot be collected by other clusters
     *  @param  clusterToCaloHitListMap the output map [parent cluster -> list of hits which belong to the main track]
     */    
    void GetExtrapolatedCaloHits(const ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList,
        const pandora::ClusterList &createdMainTrackClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    /**
     *  @brief  Perform topological checks on the collected hits to ensure no gaps are present
     *
     *  @param  clusterToCaloHitListMap the input map [parent cluster -> list of hits which belong to the main track]
     *  @param  isHigherXboundary whether the endpoint is closer to the higher x tpc boundary
     *  @param  clusterEndpointAssociation the clusterEndpointAssociation
     *
     *  @return  bool wether the checks pass
     */        
    bool AreExtrapolatedHitsGood(const ClusterToCaloHitListMap &clusterToCaloHitListMap, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation) const;

    /**
     *  @brief  Check the separation of the extremal extrapolated hits with the TPC boundary and the clusterMergePoint or, in the case of no hits, the clusterMergePoint with the TPC boundary
     *
     *  @param  clusterHitVector the extrapolated hit vector (ordered closest hit to the upstream merge point -> furthest hit)
     *  @param  isHigherXboundary whether the endpoint is closer to the higher x tpc boundary
     *  @param  clusterEndpointAssociation the clusterEndpointAssociation
     *
     *  @return  bool wether the checks pass
     */      
    bool IsExtrapolatedEndpointNearBoundary(const pandora::CaloHitVector &extrapolatedHitVector, const bool isHigherXBoundary, 
        ClusterEndpointAssociation &clusterAssociation) const;

    /**
     *  @brief  Remove the cluster association from the cluster vector so that the cluster endpoint is not considered again
     *
     *  @param  pOldConsideredCluster the cluster address of the unmodified cluster
     *  @param  pNewConsideredCluster the cluster address of the modified cluster
     *  @param  clusterVector the vector of clusters considered in future iterations of the algorithm
     *  @param  consideredClusters the list of clusters for which an endpoint has been considered
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     */       
    void ConsiderClusterAssociation(const pandora::Cluster *const pOldConsideredCluster, const pandora::Cluster *const pNewConsideredCluster,
        pandora::ClusterVector &clusterVector, pandora::ClusterList &consideredClusters, SlidingFitResultMapPair &slidingFitResultMapPair) const;    

    /**
     *  @brief  Refine the cluster endpoints and, if extrapolated hits have been found, extend to the TPC boundary and handle any remnant clusters
     *
     *  @param  clusterAssociation the cluster endpoint association
     *  @param  clusterToCaloHitListMap the map [parent cluster -> list of hits which belong to the main track]
     *  @param  pClusterList the list of all clusters
     *  @param  isHigherXboundary whether the endpoint is closer to the higher x tpc boundary
     *  @param  clusterVector the vector of clusters considered in future iterations of the algorithm
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     *  @param  consideredClusters the list of clusters for which an endpoint has been considered
     *
     *  @return  Cluster the address of the created main track cluster
     */
    const pandora::Cluster *CreateMainTrack(const ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
        const pandora::ClusterList *pClusterList, const bool isHigherXBoundary, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair,
        pandora::ClusterList &consideredClusters) const;

    unsigned int    m_maxLoopIterations;          ///< The maximum number of main loop iterations
    float           m_growingFitInitialLength;    ///< The length of hits used to initialise the extrapolated hits running fit
    float           m_growingFitSegmentLength;    ///< The length of the extrapolated hits running fit segments
    float           m_distanceToLine;             ///< The threshold hit distance of an extrapolated hit from the segment connecting line
    float           m_boundaryTolerance;          ///< The maximum allowed distance of an extremal extrapolate hit to a reference point i.e. TPC boundary
    float           m_detectorMinXEdge;           ///< The minimum x coordinate of the detector edge
    float           m_detectorMaxXEdge;           ///< The maximum x coordinate of the detector edge
    float           m_tpcMinXEdge;                ///< The minimum x coordinate of the tpc volume
    float           m_tpcMaxXEdge;                ///< The maximum x coordinate of the tpc volume
};
    
} // namespace lar_content

#endif // #ifndef LAR_TRACK_EXTENSION_REFINEMENT_ALGORITHM_H