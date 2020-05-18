/**
 *  @file   EMTrackAlgorithm.h
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
 *  @brief  EMTrackAlgorithm class
 */
class EMTrackAlgorithm : public pandora::Algorithm
{
public:

    class ClusterAssociation
    {
    public:
        ClusterAssociation();
        ClusterAssociation(const pandora::Cluster *const pInnerCluster, const pandora::Cluster *const pOuterCluster, const pandora::CartesianVector &innerMergePoint, const pandora::CartesianVector &innerMergeDirection, const pandora::CartesianVector &outerMergePoint, const pandora::CartesianVector &outerMergeDirection, const pandora::CartesianVector &connectingLineDirection);

        const pandora::Cluster* GetInnerCluster() const;
        const pandora::Cluster* GetOuterCluster() const;
        pandora::CartesianVector GetInnerMergePoint() const;
        pandora::CartesianVector GetInnerMergeDirection() const;
        pandora::CartesianVector GetOuterMergePoint() const;
        pandora::CartesianVector GetOuterMergeDirection() const;
        pandora::CartesianVector GetConnectingLineDirection() const;
    private:
        const pandora::Cluster*     m_pInnerCluster;
        const pandora::Cluster*     m_pOuterCluster;
        pandora::CartesianVector    m_innerMergePoint;
        pandora::CartesianVector    m_innerMergeDirection;
        pandora::CartesianVector    m_outerMergePoint;
        pandora::CartesianVector    m_outerMergeDirection;
        pandora::CartesianVector    m_connectingLineDirection;
    };

    /**
     *  @brief  Default constructor
     */
    EMTrackAlgorithm();

    typedef std::unordered_map<const pandora::Cluster*, const pandora::Cluster*> ClusterToAssociatedClusterMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> CaloHitToParentClusterMap;
    typedef std::unordered_map<const pandora::Cluster*, ClusterAssociation> ClusterToClusterAssociationMap;
    
 private:
    pandora::StatusCode Run();

    class SortByDistanceToLine
        {
        public:
	        SortByDistanceToLine(const pandora::CartesianVector referencePoint, const pandora::CartesianVector referenceDirection) : m_referencePoint(referencePoint), m_referenceDirection(referenceDirection.GetUnitVector()) {}
           
            bool operator() (const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs);

        private:
            const pandora::CartesianVector m_referencePoint;   ///< The point relative to which constituent hits are ordered
            const pandora::CartesianVector m_referenceDirection;
        };
    
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector);

    void InitialiseSlidingFitResultMaps(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap);

    bool FindBestClusterAssociation(const pandora::Cluster *const pCurrentCluster, const pandora::ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
        const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation);

    bool GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult, const TwoDSlidingFitResult &associatedMacroFitResult,
        pandora::CartesianVector &currentMergePosition, pandora::CartesianVector &currentMergeDirection, const bool isInner);

    bool AreClustersAssociated(const pandora::CartesianVector &currentPoint, const pandora::CartesianVector &currentDirection, const pandora::CartesianVector &testPoint, const pandora::CartesianVector &testDirection);

    void GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, pandora::CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap);

    bool IsTrackContinuous(const ClusterAssociation &clusterAssociation, pandora::CaloHitVector &extrapolatedCaloHitVector);

    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point);

    void RefineTracks(const ClusterAssociation &clusterAssociation, const pandora::CaloHitVector &extrapolatedCaloHitVector,
        const TwoDSlidingFitResultMap &microFitResultMap);

    void RefineTrack(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition, const pandora::CaloHitVector &extrapolatedCaloHitVector,
        const TwoDSlidingFitResultMap &microFitResultMap, const bool isInner);
    
    void UpdateSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap);

    void RemoveClusterFromSlidingFitResultMaps(const pandora::Cluster *const pCluster, std::vector<TwoDSlidingFitResultMap*> &slidingFitResultMapVector);
    
    void RemoveClusterFromClusterVector(const pandora::Cluster *const pCluster, pandora::ClusterVector &clusterVector);
    
    void AddHitsToCluster(const ClusterAssociation &clusterAssociation, const CaloHitToParentClusterMap &caloHitToParentClusterMap,
        const pandora::CaloHitVector &extrapolatedCaloHitVector, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
        TwoDSlidingFitResultMap &macroSlidingFitResultMap);


    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;
    CaloHitToParentClusterMap m_caloHitToParentClusterMap;
    
    float m_minCaloHits; // minimum number of calo hits for cluster to be considered
    bool m_minSeparationDistance;
    float m_maxXSeparation;
    float m_maxZSeparation;
    unsigned int m_slidingFitWindow;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* EMTrackAlgorithm::ClusterAssociation::GetInnerCluster() const
{
    return m_pInnerCluster;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* EMTrackAlgorithm::ClusterAssociation::GetOuterCluster() const
{
    return m_pOuterCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetInnerMergePoint() const
{
    return m_innerMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetInnerMergeDirection() const
{
    return m_innerMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetOuterMergePoint() const
{
    return m_outerMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetOuterMergeDirection() const
{
    return m_outerMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetConnectingLineDirection() const
{
    return m_connectingLineDirection;
}

} //namespace lar_content


#endif // #ifndef LAR_EM_TRACK_ALGORITHM_H
