/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackInEMShowerAlgorithm.h
 *
 *  @brief  Header file for the track in em shower algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_IN_EM_SHOWER_ALGORITHM_H
#define LAR_TRACK_IN_EM_SHOWER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Objects/CartesianVector.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include <unordered_map>

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
         *
         *  @return  Cluster the address of the upstream cluster
         */
        const pandora::Cluster *GetUpstreamCluster() const;

        /**
         *  @brief  Returns the downstream cluster address
         *
         *  @return  Cluster the address of the downstream cluster
         */
        const pandora::Cluster *GetDownstreamCluster() const;

        /**
         *  @brief  Returns the upstream cluster merge point
         *
         *  @return  CartesianVector the merge point of the upstream cluster
         */
        const pandora::CartesianVector &GetUpstreamMergePoint() const;

        /**
         *  @brief  Returns the upstream cluster direction at the upstream merge point
         *
         *  @return  CartesianVector the direction at the merge point of the upstream cluster
         */
        const pandora::CartesianVector &GetUpstreamMergeDirection() const;

        /**
         *  @brief  Returns the downstream cluster merge point
         *
         *  @return  CartesianVector the merge point of the downstream cluster
         */
        const pandora::CartesianVector &GetDownstreamMergePoint() const;

        /**
         *  @brief  Returns the downstream cluster direction at the downstream merge point
         *
         *  @return  CartesianVector the direction at the merge point of the downstream cluster
         */
        const pandora::CartesianVector &GetDownstreamMergeDirection() const;

        /**
         *  @brief  Returns the unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
         *
         *  @return  CartesianVector the unit displacement vector from the upstream merge point to the downstream merge point
         */
        const pandora::CartesianVector &GetConnectingLineDirection() const;

    private:
        const pandora::Cluster     *m_pUpstreamCluster;            ///< The upstream cluster of the two associated clusters
        const pandora::Cluster     *m_pDownstreamCluster;          ///< The downstream cluster of the two associated clusters
        pandora::CartesianVector    m_upstreamMergePoint;          ///< The upstream cluster point to be used in the merging process
        pandora::CartesianVector    m_upstreamMergeDirection;      ///< The upstream cluster direction at the upstream merge point (points in the direction of the downstream cluster)
        pandora::CartesianVector    m_downstreamMergePoint;        ///< The downstream cluster point to be used in the merging process
        pandora::CartesianVector    m_downstreamMergeDirection;    ///< The downstream cluster direction at the downstream merge point (points in the direction of the upstream cluster)
        pandora::CartesianVector    m_connectingLineDirection;     ///< The unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
    };

    /**
     *  @brief  Default constructor
     */
    TrackInEMShowerAlgorithm();

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
	    SortByDistanceAlongLine(const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection);

        /**
         *  @brief  Sort hits by their perpendicular distance to a line
         *
         *  @param  pLhs the address of the first hit
         *  @param  pRhs the address of the second hit
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

    typedef std::unordered_map<const pandora::Cluster*, pandora::CaloHitList> ClusterToCaloHitListMap;
    typedef std::vector<TwoDSlidingFitResultMap*> SlidingFitResultMapVector;

    /**
     *  @brief  Select clusters to be considered in algorithm
     *
     *  @param  pClusterList the input cluster list
     *  @param  clusterVector the output cluster vector
     */
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Cache the sliding fits of input clusters
     *
     *  @param  clusterVector the input cluster vector
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global cluster gradients
     */
    void InitialiseSlidingFitResultMaps(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

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
    bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
        const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation) const;

    /**
     *  @brief  Get the merging coordinate and direction for an input cluster with respect to an associated cluster
     *
     *  @param  currentMicroFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  currentMacroFitResult the global TwoDSlidingFitResult of the cluster
     *  @param  associatedMacroFitReult the global TwoDSlidingFitResult of the associated cluster
     *  @param  isUpstream whether the cluster is the upstream cluster
     *  @param  currentMergePosition the merge position of the cluster
     *  @param  currentMergeDirection the merge direction of the cluster
     *
     *  @return bool whether it was possible to find a suitable merge position
     */
    bool GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
        const TwoDSlidingFitResult &associatedMacroFitResult, const bool isUpstream, pandora::CartesianVector &currentMergePosition,
        pandora::CartesianVector &currentMergeDirection) const;

    /**
     *  @brief  Whether two clusters are assoicated to one another
     *
     *  @param  upstreamPoint the merge point of the upstream cluster
     *  @param  upstreamDirection the local direction of the upstream cluster at the merge point (in the direction of the downstream cluster)
     *  @param  downstreamPoint the merge point of the downstream cluster
     *  @param  downstreamDirection the local direction of the downstrean cluster at the merge point (in the direction of the upstream cluster)
     *
     *  @return bool whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::CartesianVector &upstreamPoint, const pandora::CartesianVector &upstreamDirection,
        const pandora::CartesianVector &downstreamPoint, const pandora::CartesianVector &downstreamDirection) const;

    /**
     *  @brief  Collect the hits that lie near the line connecting the associated clusters
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  pClusterList the list of all clusters
     *  @param  clusterToCaloHitListMap the output map [parent cluster -> list of hits which belong to the main track]
     */
    void GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList,
        ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    /**
     *  @brief  Check whether the extrapolatedCaloHitVector contains a continuous line of hits between the cluster merge points
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  clusterToCaloHitListMap the map [parent cluster -> list of hits which belong to the main track]
     *
     *  @return bool whether the calo hits are continuous
     */
    bool IsTrackContinuous(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    /**
     *  @brief  Obtain the segment boundaries of the connecting line to test whether extrapolated hits are continuous
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  trackSegmentBoundaries the output vector of segment boundaries
     */
    void GetTrackSegmentBoundaries(const ClusterAssociation &clusterAssociation, pandora::CartesianPointVector &trackSegmentBoundaries) const;

    /**
     *  @brief  Whether a calo hit falls within a specified segment of the cluster connecting line
     *
     *  @param  lowerBoundary the lower boundary of the segment
     *  @param  upperBoundary the upper boundary of the segment
     *  @param  point the position of the calo hit
     *
     *  @return bool whether the calo hit falls within segment
     */
    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    /**
     *  @brief  Merge together the split track alongside extrapolated hits
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  clusterToCaloHitListMap the map [parent cluster -> list of hits which belong to the main track]
     *  @param  pClusterList the list of all clusters
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global cluster gradients
     *  @param  clusterVector the vector of 'relevant' clusters
     */
    void CreateMainTrack(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *const pClusterList,
        TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Remove any hits in the upstream/downstream cluster that lie off of the main track axis (i.e. clustering errors)
     *
     *  @param  pCluster the input cluster
     *  @param  splitPosition the position after which hits are considered for removal
     *  @param  isUpstream whether the input cluster is the upstream cluster
     *  @param  clusterToCaloHitListMap the map [parent cluster -> list of hits which belong to the main track]
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global cluster gradients
     *  @param  clusterVector the vector of 'relevant' clusters
     *
     *  @return Cluster the address of the (possibly) modified cluster
     */
    const pandora::Cluster *RemoveOffAxisHitsFromTrack(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition, const bool isUpstream,
        const ClusterToCaloHitListMap &clusterToCaloHitListMap, pandora::ClusterList &remnantClusterList, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Remove the hits from a shower cluster that belong to the main track and add them into the main track cluster
     *
     *  @param  pShowerCluster the input shower cluster
     *  @param  pMainTrackCluster the main track cluster
     *  @param  caloHitsToMerge the list of calo hits to remove from the shower cluster
     *  @param  clusterAssociation the clusterAssociation
     *  @param  remnantClusterList the input list to store the remnant clusters
     */
    void AddHitsToMainTrack(const pandora::Cluster *const pShowerCluster, const pandora::Cluster *const pMainTrackCluster, const pandora::CaloHitList &caloHitsToMerge,
        const ClusterAssociation &clusterAssociation, pandora::ClusterList &remnantClusterList) const;

    /**
     *  @brief  Process the remnant clusters separating those that stradle the main track
     *
     *  @param  remnantClusterList the list of remnant clusters
     *  @param  pMainTrackCluster the main track cluster
     *  @param  pClusterList the list of all clusters
     */
    void ProcessRemnantClusters(const pandora::ClusterList &remnantClusterList, const pandora::Cluster *const pMainTrackCluster, const pandora::ClusterList *const pClusterList) const;

    /**
     *  @brief  Add a cluster to the nearest cluster satisfying separation distance thresholds
     *
     *  @param  pClusterToMerge the cluster to merge
     *  @param  pMainTrackCluster the main track cluster
     *  @param  pClusterList the list of all clusters
     */
    void AddToNearestCluster(const pandora::Cluster *const pClusterToMerge, const pandora::Cluster *const pClusterToEnlarge, const pandora::ClusterList *const pClusterList) const;

    /**
     *  @brief  Whether a remnant cluster is considered to be disconnected and therefore should undergo further fragmentation
     *
     *  @param  pRemnantCluster the input remnant cluster
     */
    bool IsClusterRemnantDisconnected(const pandora::Cluster *const pRemnantCluster) const;

    /**
     *  @brief  Fragment a cluster using simple hit separation logic
     *
     *  @param  pRemnantCluster the input remnant cluster to fragment
     *  @param  fragmentedClusterList the list of created clusters
     */
    void FragmentRemnantCluster(const pandora::Cluster *const pRemnantCluster, pandora::ClusterList &fragmentedClusterList) const;

    /**
     *  @brief  Update the sliding fit maps and cluster vector in preparation for a cluster deletion
     *
     *  @param  pCluster the cluster to be deleted
     *  @param  clusterVector the vector of 'relevant' clusters
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */
    void UpdateForClusterDeletion(const pandora::Cluster *const pCluster, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    /**
     *  @brief  Add only the created main track cluster to the sliding fit maps and cluster vector after the main track creation process
     *
     *  @param  pMainTrackCluster the main track cluster
     *  @param  clusterVector the vector of 'relevant' clusters
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */
    void UpdateAfterMainTrackCreation(const pandora::Cluster *const pMainTrackCluster, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    unsigned int   m_maxMainLoopIterations;                 ///< The maximum number of main loop iterations
    unsigned int   m_minCaloHits;                           ///< The threshold number of calo hits
    unsigned int   m_slidingFitWindow;                      ///< The sliding fit window used in the fits contained within the microSlidingFitResultMap
    unsigned int   m_macroSlidingFitWindow;                 ///< The sliding fit window for macro fits
    float          m_minClusterLengthSum;                   ///< The threshold cluster and associated cluster length sum
    float          m_stableRegionClusterFraction;           ///< The threshold fraction of fit contributing layers which defines the stable region
    float          m_mergePointMinCosAngleDeviation;        ///< The threshold cos opening angle between the cluster local gradient and the associated cluster global gradient used to determine merge points
    float          m_minSeparationDistance;                 ///< The threshold separation distance between associated clusters
    float          m_minDirectionDeviationCosAngle;         ///< The threshold cos opening angle of the associated cluster directions
    float          m_maxPredictedMergePointOffset;          ///< The threshold separation distance between the predicted and true cluster merge points
    float          m_distanceFromLine;                      ///< The threshold hit distance of an extrapolated hit from the cluster connecting line
    unsigned int   m_maxTrackGaps;                          ///< The maximum number of graps allowed in the extrapolated hit vector
    float          m_lineSegmentLength;                     ///< The length of a track gap
    float          m_minHitFractionForHitRemoval;           ///< The threshold fraction of hits to be removed from the cluster for hit removal to proceed
    float          m_maxDistanceFromMainTrack;              ///< The threshold distance for a hit to be added to the main track
    float          m_maxHitDistanceFromCluster;             ///< The threshold separation between a hit and cluster for the hit to be merged into the cluster
    float          m_maxHitSeparationForConnectedCluster;   ///< The maximum separation between two adjacent (in z) hits in a connected cluster
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *TrackInEMShowerAlgorithm::ClusterAssociation::GetUpstreamCluster() const
{
    return m_pUpstreamCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *TrackInEMShowerAlgorithm::ClusterAssociation::GetDownstreamCluster() const
{
    return m_pDownstreamCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TrackInEMShowerAlgorithm::ClusterAssociation::GetUpstreamMergePoint() const
{
    return m_upstreamMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TrackInEMShowerAlgorithm::ClusterAssociation::GetUpstreamMergeDirection() const
{
    return m_upstreamMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TrackInEMShowerAlgorithm::ClusterAssociation::GetDownstreamMergePoint() const
{
    return m_downstreamMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TrackInEMShowerAlgorithm::ClusterAssociation::GetDownstreamMergeDirection() const
{
    return m_downstreamMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TrackInEMShowerAlgorithm::ClusterAssociation::GetConnectingLineDirection() const
{
    return m_connectingLineDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackInEMShowerAlgorithm::SortByDistanceAlongLine::SortByDistanceAlongLine(const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection) :
    m_startPoint(startPoint),
    m_lineDirection(lineDirection.GetUnitVector())
{
}

} //namespace lar_content

#endif // #ifndef LAR_EM_TRACK_ALGORITHM_H
