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
class TrackExtensionRefinementAlgorithm :  public TrackRefinementBaseAlgorithm<ClusterEndpointAssociation>
{
public:
    
    TrackExtensionRefinementAlgorithm();

    class SortByDistanceToTPCBoundary
    {
    public:

        SortByDistanceToTPCBoundary(const float tpcXBoundary);
        bool operator() (const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    private:
        float m_tpcXBoundary;
    };
    
protected:

    typedef std::vector<ClusterEndpointAssociation> ClusterAssociaEtionVector;
    
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;
    
    virtual bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, const bool isHigherXBoundary) = 0;    

    pandora::StatusCode Run();
    
    void GetExtrapolatedCaloHits(ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    bool IsExtrapolatedEndpointNearBoundary(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const float boundaryTolerance) const;

    void CreateMainTrack(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, pandora::ClusterList &consideredClusters) const;

    void ConsiderCluster(const ClusterEndpointAssociation &clusterAssociation, pandora::ClusterVector &clusterVector) const;

    void InitialiseGeometry();


    
    bool AreExtrapolatedHitsGood(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const bool isHigherXBoundary) const;

    bool IsExtrapolatedEndpointNearBoundary(const pandora::CaloHitVector &extrapolatedHitVector, const bool isHigherXBoundary, const float boundaryTolerance, 
        ClusterEndpointAssociation &clusterAssociation) const;


    float GetAverageHitSeparation(const pandora::CaloHitVector &orderedCaloHitVector) const;

    // void UpdateAfterMainTrackModification(const pandora::Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation,
    //SlidingFitResultMapPair &slidingFitResultMapPair) const;

    float m_growingFitInitialLength;
    float m_growingFitSegmentLength;
    float m_furthestDistanceToLine;
    float m_closestDistanceToLine;
               
    float m_detectorMinXEdge;
    float m_detectorMaxXEdge;
    const pandora::LArTPC *m_pLArTPC;
    float m_tpcMinXEdge;
    float m_tpcMaxXEdge;
    
};




//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackExtensionRefinementAlgorithm::SortByDistanceToTPCBoundary::SortByDistanceToTPCBoundary(const float tpcXBoundary) :
    m_tpcXBoundary(tpcXBoundary)
{
}    
    
    
} // namespace lar_content

#endif // #ifndef LAR_TRACK_EXTENSION_REFINEMENT_ALGORITHM_H
