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
         *  @param  pInnerCluster the inner cluster of the two associated clusters
         *  @param  pOuterCluster the outer cluster of the two associated clusters
         *  @param  innerMergePoint the inner cluster point to be used in the merging process
         *  @param  innerMergeDirection the inner cluster direction at the inner merge point
         *  @param  outerMergePoint the outer cluster point to be used in the merging process
         *  @param  outerMergeDirection the outer cluster direction at the outer merge point
         */
        ClusterAssociation(const pandora::Cluster *const pInnerCluster, const pandora::Cluster *const pOuterCluster, const pandora::CartesianVector &innerMergePoint,
            const pandora::CartesianVector &innerMergeDirection, const pandora::CartesianVector &outerMergePoint, const pandora::CartesianVector &outerMergeDirection);

        /**
         *  @brief  Returns the inner cluster address
         */
        const pandora::Cluster* GetInnerCluster() const;

        /**
         *  @brief  Returns the outer cluster address
         */
        const pandora::Cluster* GetOuterCluster() const;

        /**
         *  @brief  Returns the inner cluster merge point
         */
        pandora::CartesianVector GetInnerMergePoint() const;

        /**
         *  @brief  Returns the inner cluster direction at the inner merge point
         */        
        pandora::CartesianVector GetInnerMergeDirection() const;

        /**
         *  @brief  Returns the outer cluster merge point
         */
        pandora::CartesianVector GetOuterMergePoint() const;

        /**
         *  @brief  Returns the outer cluster direction at the outer merge point
         */        
        pandora::CartesianVector GetOuterMergeDirection() const;

        /**
         *  @brief  Returns the unit vector of the line connecting the inner and outer merge points (inner -> outer)
         */           
        pandora::CartesianVector GetConnectingLineDirection() const;

        /**
         *  @brief  Sets the inner cluster address
         */
        void SetInnerCluster(const pandora::Cluster *const pCluster);

        /**
         *  @brief  Sets the outer cluster address
         */
        void SetOuterCluster(const pandora::Cluster *const pCluster);
        
    private:
        const pandora::Cluster*     m_pInnerCluster;               ///< The inner cluster of the two associated clusters         
        const pandora::Cluster*     m_pOuterCluster;               ///< The outer cluster of the two associated clusters
        pandora::CartesianVector    m_innerMergePoint;             ///< The inner cluster point to be used in the merging process
        pandora::CartesianVector    m_innerMergeDirection;         ///< The inner cluster direction at the inner merge point (points in the direction of the outer cluster)
        pandora::CartesianVector    m_outerMergePoint;             ///< The outer cluster point to be used in the merging process
        pandora::CartesianVector    m_outerMergeDirection;         ///< The outer cluster direction at the outer merge point (points in the direction of the inner cluster)
        pandora::CartesianVector    m_connectingLineDirection;     ///< The unit vector of the line connecting the inner and outer merge points (inner -> outer)
    };

    /**
     *  @brief  Default constructor
     */
    TrackInEMShowerAlgorithm();

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> CaloHitToParentClusterMap;
    
 private:
    /**
      *  @brief  SortByDistanceToLine class
      */
    class SortByDistanceToLine
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  startPoint the line start point
         *  @param  lineDirection the line direction unit vector
         */
	    SortByDistanceToLine(const pandora::CartesianVector &startPoint, const pandora::CartesianVector &lineDirection) :
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
     *  @param  isInner whether the cluster is the inner cluster
     *
     *  @return bool whether it was possible to find a suitable merge position
     */
    bool GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
        const TwoDSlidingFitResult &associatedMacroFitResult, pandora::CartesianVector &currentMergePosition, pandora::CartesianVector &currentMergeDirection,
        const bool isInner) const;

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
        CaloHitToParentClusterMap &caloHitToParentClusterMap) const;

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
     *  @brief  Remove any hits in the inner/outer cluster that are between the merge points are not found in the extrapolatedCaloHitVector
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
     *  @param  isInner whether the cluster is the inner cluster wrt the associated vector
     */
    const pandora::Cluster* RefineTrack(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition, TwoDSlidingFitResultMap &microFitResultMap,
        TwoDSlidingFitResultMap &macroFitResultMap, const bool isInner, pandora::ClusterVector &clusterVector, pandora::CaloHitVector &extrapolatedCaloHitVector) const;

    /**
     *  @brief  Remove clusters not found in an input cluster vector from sliding fit maps
     *
     *  @param  clusterVector the input vector of clusters
     *  @param  microFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to local gradients
     *  @param  macroFitResultMap the mapping [cluster -> TwoDSlidingFitResult] where fits correspond to global gradients
     */
    void UpdateSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;

    void UpdateForClusterCreation(const pandora::Cluster *&pCluster, PandoraContentApi::Cluster::Parameters &clusterParameters, pandora::ClusterVector &clusterVector) const;
    
    /**
     *  @brief  Remove a cluster from the sliding fit result maps in an input vector
     *
     *  @param  pCluster the input cluster
     *  @param  slidingFitResultMapVector input vector of sliding fit result maps
     */
    void RemoveClusterFromSlidingFitResultMaps(const pandora::Cluster *const pCluster, std::vector<TwoDSlidingFitResultMap*> &slidingFitResultMapVector) const;

    /**
     *  @brief  Remove a cluster from an cluster vector
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
    void AddHitsToCluster(const ClusterAssociation &clusterAssociation, const CaloHitToParentClusterMap &caloHitToParentClusterMap,
        const pandora::CaloHitVector &extrapolatedCaloHitVector, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap) const;
    
    unsigned int m_minCaloHits;                                ///< The threshold number of calo hits 
    float m_minSeparationDistance;                             ///< The threshold separation distance between associated clusters 
    float m_maxXSeparation;                                    ///< The threshold separation distance in x between extrapolated and true cluster merge points 
    float m_maxZSeparation;                                    ///< The threshold separation distance in z between extrapolated and true cluster merge points 
    unsigned int m_slidingFitWindow;                           ///< The sliding fit window used in the fits contained within the microSlidingFitResultMap
    float m_mergePointMinCosAngleDeviation;                    ///< The threshold cos opening angle between the cluster local gradient and the associated
                                                               ///  cluster global gradient used to determine merge points
    float m_minClusterLengthSum;                               ///< The threshold cluster and associated cluster length sum
    float m_minDirectionDeviationCosAngle;                     ///< The threshold cos opening angle of the associated cluster directions
    float m_distanceFromLine;                                  ///< The threshold hit distance of an extrapolated hit from the cluster connecting line 
    unsigned int m_maxTrackGaps;                               ///< The maximum number of graps allowed in the extrapolated hit vector
    float m_lineSegmentLength;                                 ///< The length of a track gap

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* TrackInEMShowerAlgorithm::ClusterAssociation::GetInnerCluster() const
{
    return m_pInnerCluster;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* TrackInEMShowerAlgorithm::ClusterAssociation::GetOuterCluster() const
{
    return m_pOuterCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetInnerMergePoint() const
{
    return m_innerMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetInnerMergeDirection() const
{
    return m_innerMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetOuterMergePoint() const
{
    return m_outerMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetOuterMergeDirection() const
{
    return m_outerMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector TrackInEMShowerAlgorithm::ClusterAssociation::GetConnectingLineDirection() const
{
    return m_connectingLineDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackInEMShowerAlgorithm::ClusterAssociation::SetInnerCluster(const pandora::Cluster *const pCluster)
{
    // ISOBEL DO I ADD 'IF CLUSTER IS NOT NULL'
    m_pInnerCluster = pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackInEMShowerAlgorithm::ClusterAssociation::SetOuterCluster(const pandora::Cluster *const pCluster)
{
    // ISOBEL DO I ADD 'IF CLUSTER IS NOT NULL'    
    m_pOuterCluster = pCluster;
}
    
} //namespace lar_content

#endif // #ifndef LAR_EM_TRACK_ALGORITHM_H
