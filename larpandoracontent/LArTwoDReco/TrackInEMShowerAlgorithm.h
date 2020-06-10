/**
 *  @file   TrackInEMShowerAlgorithm.h
 *
 *  @brief  Header file for the em track algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_EM_TRACK_ALGORITHM_H
#define LAR_EM_TRACK_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  TrackInEMShowerAlgorithm class
 */
class TrackInEMShowerAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  ClusterAssociation class
     */
    class ClusterAssociation
    {
    public:
        /**
         *  @brief  Default constructor
         */
        ClusterAssociation();

        /**
         *  @brief  Constructor
         *
         *  @param  pUpstreamCluster the upstream cluster of the two associated clusters
         *  @param  pDownstreamCluster the downstream cluster of the two associated clusters
         *  @param  upstreamMergePoint the upstream cluster point to be used in the merging process
         *  @param  upstreamMergeDirection the upstream cluster direction at the upstream merge point
         *  @param  downstreamMergePoint the downstream cluster point to be used in the merging process
         *  @param  downstreamMergeDirection the downstream cluster direction at the downstream merge point
         */
        ClusterAssociation(const pandora::Cluster *const pUpstreamCluster, const pandora::Cluster *const pDownstreamCluster, const pandora::CartesianVector &upstreamMergePoint,
            const pandora::CartesianVector &upstreamMergeDirection, const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection);

        /**
         *  @brief  Returns the upstream cluster address
         */
        const pandora::Cluster* GetUpstreamCluster() const;

        /**
         *  @brief  Returns the downstream cluster address
         */
        const pandora::Cluster* GetDownstreamCluster() const;

        /**
         *  @brief  Returns the upstream cluster merge point
         */
        pandora::CartesianVector GetUpstreamMergePoint() const;

        /**
         *  @brief  Returns the upstream cluster direction at the upstream merge point
         */        
        pandora::CartesianVector GetUpstreamMergeDirection() const;

        /**
         *  @brief  Returns the downstream cluster merge point
         */
        pandora::CartesianVector GetDownstreamMergePoint() const;

        /**
         *  @brief  Returns the downstream cluster direction at the downstream merge point
         */        
        pandora::CartesianVector GetDownstreamMergeDirection() const;

        /**
         *  @brief  Returns the unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
         */           
        pandora::CartesianVector GetConnectingLineDirection() const;

        /**
         *  @brief  Sets the upstream cluster address
         */
        void SetUpstreamCluster(const pandora::Cluster *const pCluster);

        /**
         *  @brief  Sets the downstream cluster address
         */
        void SetDownstreamCluster(const pandora::Cluster *const pCluster);
        
    private:
        const pandora::Cluster*     m_pUpstreamCluster;               ///< The upstream cluster of the two associated clusters         
        const pandora::Cluster*     m_pDownstreamCluster;               ///< The downstream cluster of the two associated clusters
        pandora::CartesianVector    m_upstreamMergePoint;             ///< The upstream cluster point to be used in the merging process
        pandora::CartesianVector    m_upstreamMergeDirection;         ///< The upstream cluster direction at the upstream merge point (points in the direction of the downstream cluster)
        pandora::CartesianVector    m_downstreamMergePoint;             ///< The downstream cluster point to be used in the merging process
        pandora::CartesianVector    m_downstreamMergeDirection;         ///< The downstream cluster direction at the downstream merge point (points in the direction of the upstream cluster)
        pandora::CartesianVector    m_connectingLineDirection;     ///< The unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
    };

    /**
     *  @brief  Default constructor
     */
    TrackInEMShowerAlgorithm();

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> CaloHitToParentClusterMap;
    typedef std::unordered_map<const pandora::Cluster*, pandora::CaloHitList> ClusterToCaloHitListMap;
    
 private:
    /**
      *  @brief  SortByDistanceAlongLine class
      */
    class SortByDistanceAlongLine
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  startPoint the line start point
         *  @param  lineDirection the line direction unit vector
         */
	    SortByDistanceAlongLine(const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection) :
            m_startPoint(startPoint), m_lineDirection(lineDirection.GetUnitVector()) {}

        /**
         *  @brief  Sort hits by their perpendicular distance to a line
             *
             *  @param  lhs first hit
             *  @param  rhs second hit
             *
             *  @return  bool
             */
        bool operator() (const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs);

    private:
        const pandora::CartesianVector m_startPoint;        ///< The line start point
        const pandora::CartesianVector m_lineDirection;     ///< The line end point
    };

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    /**
     *  @brief  Select clusters to be considered in algorithm
     *
     *  @param  pClusterList the input cluster list
     *  @param  clusterVector the output cluster vector
     *
     */
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Cache the sliding fits of input clusters
     *
     *  @param  clusterVector the input cluster vector
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global cluster gradients
     *
     */
    void InitialiseSlidingFitResultMaps(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    /**
     *  @brief  Find the cluster association with the longest cluster length sum
     *
     *  @param  clusterVector the input cluster vector
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global cluster gradients
     *  @param  clusterAssociation the output ClusterAssociation
     *
     *  @return bool whether a cluster association was found
     */
    bool FindBestClusterAssociation(pandora::ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
        const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation) const;

    /**
     *  @brief  Get the merging coordinate and direction for an input cluster with respect to an associated cluster
     *
     *  @param  currentMicroFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  currentMacroFitResult the global TwoDSlidingFitResult of the cluster
     *  @param  associatedMacroFitReult the global TwoDSlidingFitResult of the associated cluster
     *  @param  currentMergePosition the merge position of the cluster
     *  @param  currentMergeDirection the merge direction of the cluster
     *  @param  isUpstream whether the cluster is the upstream cluster
     *
     *  @return bool whether it was possible to find a suitable merge position
     */
    bool GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
        const TwoDSlidingFitResult &associatedMacroFitResult, pandora::CartesianVector &currentMergePosition, pandora::CartesianVector &currentMergeDirection,
        const bool isUpstream) const;

    /**
     *  @brief  Whether two clusters are assoicated to one another
     *
     *  @param  currentPoint the merge point of the first cluster
     *  @param  currentDirection the local direction of the first cluster at the merge point (in the direction of the second cluster)
     *  @param  testPoint the merge point of the second cluster
     *  @param  testDirection the local direction of the second cluster at the merge point (in the direction of the first cluster)
     *
     *  @return bool whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::CartesianVector &currentPoint, const pandora::CartesianVector &currentDirection, const pandora::CartesianVector &testPoint,
        const pandora::CartesianVector &testDirection) const;

    /**
     *  @brief  Collect the hits that lie near the line connecting the associated clusters
     *
     *  @param  clusterAssociation the clusterAssociation 
     *  @param  pClusterList the list of all clusters
     *  @param  extrapolatedCaloHitVector the output vector of collected hits
     *  @param  caloHitToParentClusterMap the output map [calo hit -> parent cluster] 
     */
    void GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, pandora::CaloHitVector &extrapolatedCaloHitVector,
         ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    /**
     *  @brief  Check whether the extrapolatedCaloHitVector contains a continuous line of hits between the cluster merge points
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  extrapolatedCaloHitVector the input vector of calo hits
     *
     *  @return bool whether the calo hits are continuous
     */
    bool IsTrackContinuous(const ClusterAssociation &clusterAssociation, pandora::CaloHitVector &extrapolatedCaloHitVector) const;

    /**
     *  @brief  Whether a calo hit falls within a specified segment of the cluster connecting line
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  extrapolatedCaloHitVector the input vector of calo hits
     *
     *  @return bool whether the calo hits are continuous
     */
    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    /**
     *  @brief  Remove any hits in the upstream/downstream cluster that are between the merge points are not found in the extrapolatedCaloHitVector
     *
     *  @param clusterAssociation the clusterAssociation
     *  @param  extrapolatedCaloHitVector the input vector of calo hits
     *  @param  microFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     */
    void RefineTracks(ClusterAssociation &clusterAssociation, TwoDSlidingFitResultMap &microFitResultMap,
        TwoDSlidingFitResultMap &macroFitResultMap, pandora::ClusterVector &clusterVector, pandora::CaloHitVector &extrapolatedCaloHitVector) const;

    /**
     *  @brief  Remove cluster hits that fall after the split position
     *
     *  @param  pCluster the input cluster
     *  @param  splitPosition the position the cluster after which hits are removed
     *  @param  extrapolatedCaloHitVector the input vector of calo hits
     *  @param  microFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     *  @param  isUpstream whether the cluster is the upstream cluster wrt the associated vector
     */
    const pandora::Cluster* RefineTrack(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition, TwoDSlidingFitResultMap &microFitResultMap,
        TwoDSlidingFitResultMap &macroFitResultMap, const bool isUpstream, pandora::ClusterVector &clusterVector, pandora::CaloHitVector &extrapolatedCaloHitVector) const;

    /**
     *  @brief  Remove clusters not found in an input cluster vector from sliding fit maps
     *
     *  @param  clusterVector the input vector of clusters
     *  @param  microFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */
    void UpdateSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    /**
     *  @brief  Cluster calo hits into groups motivated by their spatial separation
     *
     *  @param  caloHitList the input calo hit list
     *  @param  clusterVector the vector of 'relevant' clusters
     */
    bool CreateClusters(const pandora::CaloHitList &caloHitList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Update the sliding fit maps and cluster vector after a cluster modification
     *
     *  @param  pCluster the modified cluster
     *  @param  clusterVector the vector of 'relevant' clusters
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */    
    void UpdateAfterClusterModification(const pandora::Cluster *const pCluster, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    /**
     *  @brief  Update the sliding fit maps and cluster vector in preparation for a cluster deletion
     *
     *  @param  pCluster the modified cluster
     *  @param  clusterVector the vector of 'relevant' clusters
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */    
    void UpdateForClusterDeletion(const pandora::Cluster *const pCluster, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;
    
    /**
     *  @brief  Remove a cluster from the sliding fit result maps
     *
     *  @param  pCluster the input cluster
     *  @param  slidingFitResultMapVector input vector of sliding fit result maps
     */
    void RemoveClusterFromSlidingFitResultMaps(const pandora::Cluster *const pCluster, std::vector<TwoDSlidingFitResultMap*> &slidingFitResultMapVector) const;

    /**
     *  @brief  Remove a cluster from a cluster vector
     *
     *  @param  pCluster the input cluster
     *  @param  clusterVector the input cluster vector
     */
    void RemoveClusterFromClusterVector(const pandora::Cluster *const pCluster, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Merge the clusters in the clusterAssociation together alongside the hits within the extrapolatedCaloHitVector. In this process 'deleted' clusters are removed from the 
     *          cached sliding fit result maps and cluster vector but left in the caloHitToParentClusterMap
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  caloHitToParentClusterMap the map [calo hit -> parent cluster] 
     *  @param  extrapolatedCaloHitVector the vector of extrapolated calo hits to be merged in to the resulting cluster
     *  @param  clusterVector the cluster vector
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */
    void AddHitsToCluster(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
        pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    void FragmentCluster(const pandora::Cluster *const pCluster, const pandora::Cluster *const pClusterToEnlarge, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
        pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;
    
    unsigned int m_minCaloHits;                                ///< The threshold number of calo hits 
    float m_minSeparationDistance;                             ///< The threshold separation distance between associated clusters 
    float m_maxPredictedMergePointOffset;                      ///< The threshold separation distance between the predicted and true cluster merge points 
    unsigned int m_slidingFitWindow;                           ///< The sliding fit window used in the fits contained within the microSlidingFitResultMap
    float m_mergePointMinCosAngleDeviation;                    ///< The threshold cos opening angle between the cluster local gradient and the associated
                                                               ///  cluster global gradient used to determine merge points
    float m_minClusterLengthSum;                               ///< The threshold cluster and associated cluster length sum
    float m_minDirectionDeviationCosAngle;                     ///< The threshold cos opening angle of the associated cluster directions
    float m_distanceFromLine;                                  ///< The threshold hit distance of an extrapolated hit from the cluster connecting line 
    unsigned int m_maxTrackGaps;                               ///< The maximum number of graps allowed in the extrapolated hit vector
    float m_lineSegmentLength;                                 ///< The length of a track gap
    float m_maxHitDistanceFromCluster;                         ///< The threshold separation between a hit and cluster for the hit to be merged into the cluster
    unsigned int m_globalSlidingFitWindow;
    
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* TrackInEMShowerAlgorithm::ClusterAssociation::GetUpstreamCluster() const
{
    return m_pUpstreamCluster;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* TrackInEMShowerAlgorithm::ClusterAssociation::GetDownstreamCluster() const
{
    return m_pDownstreamCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetUpstreamMergePoint() const
{
    return m_upstreamMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetUpstreamMergeDirection() const
{
    return m_upstreamMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetDownstreamMergePoint() const
{
    return m_downstreamMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetDownstreamMergeDirection() const
{
    return m_downstreamMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetConnectingLineDirection() const
{
    return m_connectingLineDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackInEMShowerAlgorithm::ClusterAssociation::SetUpstreamCluster(const pandora::Cluster *const pCluster)
{
    // ISOBEL DO I ADD 'IF CLUSTER IS NOT NULL'
    m_pUpstreamCluster = pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackInEMShowerAlgorithm::ClusterAssociation::SetDownstreamCluster(const pandora::Cluster *const pCluster)
{
    // ISOBEL DO I ADD 'IF CLUSTER IS NOT NULL'    
    m_pDownstreamCluster = pCluster;
}
    
} //namespace lar_content

#endif // #ifndef LAR_EM_TRACK_ALGORITHM_H
