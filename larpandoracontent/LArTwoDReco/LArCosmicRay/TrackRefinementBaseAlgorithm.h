/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/TrackRefinementBaseAlgorithm.h
 *
 *  @brief  Header file for the track refinement base class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_REFINEMENT_BASE_ALGORITHM_H
#define LAR_TRACK_REFINEMENT_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArTwoDReco/LArCosmicRay/ClusterAssociation.h"

namespace lar_content
{
/**
 *  @brief TrackRefinementBaseAlgorithm class
 */
class TrackRefinementBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackRefinementBaseAlgorithm();

protected:
    typedef std::pair<TwoDSlidingFitResultMap *, TwoDSlidingFitResultMap *> SlidingFitResultMapPair;
    typedef std::unordered_map<const pandora::Cluster *, pandora::CaloHitList> ClusterToCaloHitListMap;

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
         *  @param  hitWidthMode whether to consider hit widths or not
         */
        SortByDistanceAlongLine(const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection, const bool hitWidthMode);

        /**
         *  @brief  Sort hits by their projected distance along a line from a start point
         *
         *  @param  pLhs the address of the first hit
         *  @param  pRhs the address of the second hit
         *
         *  @return  whether lhs hit has a smaller projected distance along the line than rhs hit
         */
        bool operator()(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs) const;

    private:
        pandora::CartesianVector m_startPoint;    ///< The line start point
        pandora::CartesianVector m_lineDirection; ///< The line end point
        bool m_hitWidthMode;                      ///< Wether to consider hit widths or not
    };

    virtual pandora::StatusCode Run() = 0;
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;

    /**
     *  @brief  Fill the cluster vector and sliding fit maps with clusters that are determined to be track-like
     *
     *  @param  pClusterList the list of input clusters
     *  @param  sortFunction the sort class or function used to sort the clusterVector
     *  @param  clusterVector the input vector to store clusters considered within the algorithm
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps
     */
    template <typename T>
    void InitialiseContainers(const pandora::ClusterList *pClusterList, const T sortFunction, pandora::ClusterVector &clusterVector,
        SlidingFitResultMapPair &slidingFitResultMapPair) const;

    /**
     *  @brief  Get the merging coordinate and direction for an input cluster with respect to an associated cluster
     *
     *  @param  clusterMicroFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  clusterMacroFitResult the global TwoDSlidingFitResult of the cluster
     *  @param  associatedMacroFitReult the global TwoDSlidingFitResult of the associated cluster
     *  @param  isEndUpstream whether the sought cluster merge point is the upstream
     *  @param  clusterMergePosition the merge position of the cluster
     *  @param  clusterMergeDirection the merge direction of the cluster
     *
     *  @return  whether it was possible to find a suitable merge position
     */
    bool GetClusterMergingCoordinates(const TwoDSlidingFitResult &clusterMicroFitResult, const TwoDSlidingFitResult &clusterMacroFitResult,
        const TwoDSlidingFitResult &associatedMacroFitResult, const bool isEndUpstream, pandora::CartesianVector &clusterMergePosition,
        pandora::CartesianVector &clusterMergeDirection) const;

    /**
     *  @brief  Find the unprotected hits that are contained within a defined box with the option to apply a cut on the distance to the connecting line
     *
     *  @param  firstCorner the position of one corner
     *  @param  secondCorner the position of the opposite corner
     *  @param  pClusterList the list of all clusters
     *  @param  clusterToCaloHitListMap the output map [parent cluster -> list of hits which belong to the main track]
     *  @param  unavailableProtectedClusters the list of clusters whose hits are protected
     *  @param  distanceToLine the maximum perpendicular distance of a collected hit from the connecting line
     */
    void GetHitsInBoundingBox(const pandora::CartesianVector &firstCorner, const pandora::CartesianVector &secondCorner,
        const pandora::ClusterList *const pClusterList, ClusterToCaloHitListMap &clusterToCaloHitListMap,
        const pandora::ClusterList &unavailableProtectedClusters = pandora::ClusterList(), const float distanceToLine = -1.f) const;

    /**
     *  @brief  check whether a hit is contained within a defined square region
     *
     *  @param  minX the minimum x coordinate of the square region
     *  @param  maxX the maximum x coordinate of the square region
     *  @param  minZ the minimum z coordinate of the square region
     *  @param  maxZ the maximum z coordinate of the square region
     *  @param  hitPosition the position of the hit
     *
     *  @return  whether the hit is contained within the square region
     */
    bool IsInBoundingBox(const float minX, const float maxX, const float minZ, const float maxZ, const pandora::CartesianVector &hitPosition) const;

    /**
     *  @brief  Check whether a hit is close to a line
     *
     *  @param  hitPosition the position of the hit
     *  @param  lineStart the start point of the line (can actually be any point on the line)
     *  @param  lineDirection the unit vector of the line direction
     *  @param  distanceToLine the definition of 'close'
     *
     *  @return  whether the hit is close to the line
     */
    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart,
        const pandora::CartesianVector &lineDirection, const float distanceToLine) const;

    /**
     *  @brief  Perform topological checks on the collected hits to ensure no gaps are present
     *
     *  @param  clusterToCaloHitListMap the input map [parent cluster -> list of hits which belong to the main track]
     *  @param  clusterAssociation the clusterAssociation
     *
     *  @return  whether the checks pass
     */
    bool AreExtrapolatedHitsGood(const ClusterToCaloHitListMap &clusterToCaloHitListMap, ClusterAssociation &clusterAssociation) const;

    /**
     *  @brief  Check the separation of the extremal extrapolated hits with the expected endpoints or, in the case of no hits, the clusterMergePoint themselves
     *
     *  @param  extrapolatedHitVector the extrapolated hit vector (ordered closest hit to the upstream merge point -> furthest hit)
     *  @param  clusterAssociation the clusterAssociation
     *
     *  @return  whether the check passes
     */
    virtual bool AreExtrapolatedHitsNearBoundaries(const pandora::CaloHitVector &extrapolatedHitVector, ClusterAssociation &clusterAssociation) const = 0;

    /**
     *  @brief  Check whether a hit is close to a boundary point
     *
     *  @param  pCaloHit the input calo hit
     *  @param  boundaryPosition2D the position of the 2D boundary point
     *  @param  boundaryTolerance the definition of close
     *
     *  @return  whether the check passes
     */
    bool IsNearBoundary(const pandora::CaloHit *const pCaloHit, const pandora::CartesianVector &boundaryPosition2D, const float boundaryTolerance) const;

    /**
     *  @brief  Check whether the extrapolatedCaloHitVector contains a continuous line of hits between the cluster merge points
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  extrapolatedCaloHitVector the vector of extrapolated calo hits
     *
     *  @return  whether the calo hits are continuous
     */
    bool IsTrackContinuous(const ClusterAssociation &clusterAssociation, const pandora::CaloHitVector &extrapolatedCaloHitVector) const;

    /**
     *  @brief  Obtain the segment boundaries of the connecting line to test whether extrapolated hits are continuous
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  trackSegmentBoundaries the output vector of segment boundaries
     */
    void GetTrackSegmentBoundaries(const ClusterAssociation &clusterAssociation, pandora::CartesianPointVector &trackSegmentBoundaries) const;

    /**
     *  @brief  Move an input position to the higher line gap edge if it lies within a gap
     *
     *  @param  mergeDirection the direction of the track
     *  @param  trackPoint the input position
     */
    void RepositionIfInGap(const pandora::CartesianVector &mergeDirection, pandora::CartesianVector &trackPoint) const;

    /**
     *  @brief  Calculate the track length between two points that lies in gaps
     *
     *  @param  upstreamPoint the upstream point
     *  @param  downstreamPoint the downstream point
     *  @param  connectingLine the track direction
     *  @param  consideredGaps the list of gaps to ignore
     */
    float DistanceInGap(const pandora::CartesianVector &upstreamPoint, const pandora::CartesianVector &downstreamPoint,
        const pandora::CartesianVector &connectingLine, pandora::DetectorGapList &consideredGaps) const;

    /**
     *  @brief  Whether a position falls within a specified segment of the cluster connecting line
     *
     *  @param  lowerBoundary the lower boundary of the segment
     *  @param  upperBoundary the upper boundary of the segment
     *  @param  point the position
     *
     *  @return  whether the position falls within segment
     */
    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary,
        const pandora::CartesianVector &point) const;

    /**
     *  @brief  Remove any hits in the upstream/downstream cluster that lie off of the main track axis (i.e. clustering errors)
     *
     *  @param  pCluster the input cluster
     *  @param  splitPosition the position after which hits are considered for removal
     *  @param  isEndUpstream whether the upstream end is to be refined
     *  @param  clusterToCaloHitListMap the map [parent cluster -> list of hits which belong to the main track]
     *  @param  remnantClusterList the input list to store the remnant clusters
     *  @param  microSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroSlidingFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global cluster gradients
     *
     *  @return  the address of the (possibly) modified cluster
     */
    const pandora::Cluster *RemoveOffAxisHitsFromTrack(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition,
        const bool isEndUpstream, const ClusterToCaloHitListMap &clusterToCaloHitListMap, pandora::ClusterList &remnantClusterList,
        TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    /**
     *  @brief  Remove the hits from a shower cluster that belong to the main track and add them into the main track cluster
     *
     *  @param  pMainTrackCluster the main track cluster
     *  @param  pShowerCluster the input shower cluster
     *  @param  caloHitsToMerge the list of calo hits to remove from the shower cluster
     *  @param  clusterAssociation the clusterAssociation
     *  @param  remnantClusterList the input list to store the remnant clusters
     */
    void AddHitsToMainTrack(const pandora::Cluster *const pMainTrackCluster, const pandora::Cluster *const pShowerTrackCluster,
        const pandora::CaloHitList &caloHitsToMerge, const ClusterAssociation &clusterAssociation, pandora::ClusterList &remnantClusterList) const;

    /**
     *  @brief  Process the remnant clusters separating those that stradle the main track
     *
     *  @param  remnantClusterList the list of remnant clusters to process
     *  @param  pMainTrackCluster the main track cluster
     *  @param  pClusterList the list of all clusters
     *  @param  createdClusters the input list to store the final remnant clusters
     */
    void ProcessRemnantClusters(const pandora::ClusterList &remnantClusterList, const pandora::Cluster *const pMainTrackCluster,
        const pandora::ClusterList *const pClusterList, pandora::ClusterList &createdClusters) const;

    /**
     *  @brief  Add a cluster to the nearest cluster satisfying separation distance thresholds
     *
     *  @param  pClusterToMerge the cluster to merge
     *  @param  pMainTrackCluster the main track cluster
     *  @param  pClusterList the list of all clusters
     *
     *  @return  whether cluster was added to a nearby cluster
     */
    bool AddToNearestCluster(const pandora::Cluster *const pClusterToMerge, const pandora::Cluster *const pMainTrackCluster,
        const pandora::ClusterList *const pClusterList) const;

    /**
     *  @brief  Whether a remnant cluster is considered to be disconnected and therefore should undergo further fragmentation
     *
     *  @param  pRemnantCluster the input remnant cluster
     *
     *  @return  whether the remnant cluster is disconnected
     */
    bool IsClusterRemnantDisconnected(const pandora::Cluster *const pRemnantCluster) const;

    /**
     *  @brief  Fragment a cluster using simple hit separation logic
     *
     *  @param  pRemnantCluster the input remnant cluster to fragment
     *  @param  fragmentedClusterList the input list to store the final remnant clusters
     */
    void FragmentRemnantCluster(const pandora::Cluster *const pRemnantCluster, pandora::ClusterList &fragmentedClusterList) const;

    /**
     *  @brief  Remove deleted clusters from the cluster vector and sliding fit maps and add in created clusters that are determined to be track-like
     *
     *  @param  clustersToAdd the list of clusters to add to the containers
     *  @param  clustersToDelete the list of clusters to remove from the containers
     *  @param  sortFunction the sort class or function used to sort the clusterVector
     *  @param  clusterVector the vector to store clusters considered within the algorithm
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps
     */
    template <typename T>
    void UpdateContainers(const pandora::ClusterList &clustersToAdd, const pandora::ClusterList &clustersToDelete, const T sortFunction,
        pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    /**
     *  @brief  Remove a cluster from the cluster vector and sliding fit maps
     *
     *  @param  pClustertoRemove the clusters to remove from the containers
     *  @param  clusterVector the vector to store clusters considered within the algorithm
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps
     */
    void RemoveClusterFromContainers(const pandora::Cluster *const pClustertoRemove, pandora::ClusterVector &clusterVector,
        SlidingFitResultMapPair &slidingFitResultMapPair) const;

    float m_minClusterLength;               ///< The minimum length of a considered cluster
    unsigned int m_microSlidingFitWindow;   ///< The sliding fit window used in the fits contained within the microSlidingFitResultMap
    unsigned int m_macroSlidingFitWindow;   ///< The sliding fit window used in the fits contained within the macroSlidingFitResultMap
    float m_stableRegionClusterFraction;    ///< The threshold fraction of fit contributing layers which defines the stable region
    float m_mergePointMinCosAngleDeviation; ///< The threshold cos opening angle between the cluster local gradient and the associated cluster global gradient used to determine merge points
    float m_minHitFractionForHitRemoval;    ///< The threshold fraction of hits to be removed from the cluster for hit removal to proceed
    float m_maxDistanceFromMainTrack;       ///< The threshold distance for a hit to be added to the main track
    float m_maxHitDistanceFromCluster; ///< The threshold separation between a hit and cluster for the hit to be merged into the cluster
    float m_maxHitSeparationForConnectedCluster; ///< The maximum separation between two adjacent (in z) hits in a connected cluster
    unsigned int m_maxTrackGaps;                 ///< The maximum number of graps allowed in the extrapolated hit vector
    float m_lineSegmentLength;                   ///< The length of a track gap
    bool m_hitWidthMode;                         ///< Whether to consider the width of hits
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackRefinementBaseAlgorithm::SortByDistanceAlongLine::SortByDistanceAlongLine(
    const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection, const bool hitWidthMode) :
    m_startPoint(startPoint), m_lineDirection(lineDirection.GetUnitVector()), m_hitWidthMode(hitWidthMode)
{
}

} // namespace lar_content

#endif // #ifndef TRACK_REFINEMENT_BASE_ALGORITHM_H
