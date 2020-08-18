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

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/ClusterAssociation.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{
/**
 *  @brief TrackRefinementBaseAlgorithm class
 */
template<typename T>    
class TrackRefinementBaseAlgorithm : public pandora::Algorithm
{
public:
    
    /**
     *  @brief  Default constructor
     */
    TrackRefinementBaseAlgorithm();

protected:

    typedef std::pair<TwoDSlidingFitResultMap*, TwoDSlidingFitResultMap*> SlidingFitResultMapPair;
    typedef std::unordered_map<const pandora::Cluster*, pandora::CaloHitList> ClusterToCaloHitListMap;

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
    
    virtual pandora::StatusCode Run() = 0;
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;
    
    virtual void GetExtrapolatedCaloHits(T &clusterAssociation, const pandora::ClusterList *const pClusterList, const pandora::ClusterList &consideredClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const = 0;

    virtual void ConsiderClusterAssociation(const T &clusterAssociation, pandora::ClusterVector &clusterVector, pandora::ClusterList &consideredClusters,
        SlidingFitResultMapPair &slidingFitResultMapPair) const = 0;
    
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


    bool IsTrackContinuous(const ClusterAssociation &clusterAssociation, const pandora::CaloHitVector &extrapolatedCaloHitVector, const unsigned int maxTrackGaps,
        const float lineSegmentLength) const;

    /**
     *  @brief  Obtain the segment boundaries of the connecting line to test whether extrapolated hits are continuous
     *
     *  @param  clusterAssociation the clusterAssociation
     *  @param  trackSegmentBoundaries the output vector of segment boundaries
     */
    void GetTrackSegmentBoundaries(const ClusterAssociation &clusterAssociation, pandora::CartesianPointVector &trackSegmentBoundaries, const float lineSegmentLength) const;

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
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

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
    void ProcessRemnantClusters(const pandora::ClusterList &remnantClusterList, const pandora::Cluster *const pMainTrackCluster, const pandora::ClusterList *const pClusterList, pandora::ClusterList &createdClusters) const;

    /**
     *  @brief  Add a cluster to the nearest cluster satisfying separation distance thresholds
     *
     *  @param  pClusterToMerge the cluster to merge
     *  @param  pMainTrackCluster the main track cluster
     *  @param  pClusterList the list of all clusters
     */
    bool AddToNearestCluster(const pandora::Cluster *const pClusterToMerge, const pandora::Cluster *const pClusterToEnlarge, const pandora::ClusterList *const pClusterList) const;

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


    void InitialiseContainers(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    void UpdateContainers(const pandora::ClusterList &clustersToAdd, const pandora::ClusterList &clustersToDelete, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    // Think I'll just put this back in the other function?
    void RemoveClusterFromContainers(const pandora::Cluster *const pClustertoRemove, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    float GetAverageDeviationFromLine(const pandora::Cluster *const pCluster, const pandora::CartesianVector &line, const pandora::CartesianVector &startPoint) const;


//------------------------------------------------------------------------------------------------------------------------------------------    

    unsigned int m_minCaloHits;
    float m_maxCurviness;
    unsigned int m_microSlidingFitWindow;
    unsigned int m_macroSlidingFitWindow;
    float m_stableRegionClusterFraction;
    float m_mergePointMinCosAngleDeviation;
    float m_distanceFromLine;
    float          m_minHitFractionForHitRemoval;           ///< The threshold fraction of hits to be removed from the cluster for hit removal to proceed
    float m_maxDistanceFromMainTrack;
    float m_maxHitDistanceFromCluster;
    float m_maxHitSeparationForConnectedCluster;
    unsigned int   m_maxTrackGaps;                          ///< The maximum number of graps allowed in the extrapolated hit vector
    float          m_lineSegmentLength;                     ///< The length of a track gap
    

    
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
inline TrackRefinementBaseAlgorithm<T>::SortByDistanceAlongLine::SortByDistanceAlongLine(const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection) :
    m_startPoint(startPoint),
    m_lineDirection(lineDirection.GetUnitVector())
{
}
    
} // namespace lar_content

#endif // #ifndef TRACK_REFINEMENT_BASE_ALGORITHM_H
