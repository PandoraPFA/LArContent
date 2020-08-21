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
     *  @brief  Initialise the detector and TPC edge member parameters taking into account the gap in the case of a CPA
     */       
    void InitialiseGeometry();
    
    virtual bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        const pandora::ClusterList *const pClusterList, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation) = 0;    


    
    void GetExtrapolatedCaloHits(ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, const pandora::ClusterList &consideredClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    const pandora::Cluster *CreateMainTrack(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, pandora::ClusterList &consideredClusters, bool isHigherXBoundary) const;

    void ConsiderClusterAssociation(const pandora::Cluster *const pOldConsideredCluster, const pandora::Cluster *const pNewConsideredCluster, pandora::ClusterVector &clusterVector, pandora::ClusterList &consideredClusters, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    


    
    bool AreExtrapolatedHitsGood(ClusterEndpointAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const bool isHigherXBoundary) const;

    bool IsExtrapolatedEndpointNearBoundary(const pandora::CaloHitVector &extrapolatedHitVector, const bool isHigherXBoundary, 
        ClusterEndpointAssociation &clusterAssociation) const;


    float m_growingFitInitialLength;
    float m_growingFitSegmentLength;
    float m_furthestDistanceToLine;
    float m_closestDistanceToLine;
               
    float m_detectorMinXEdge;
    float m_detectorMaxXEdge;
    const pandora::LArTPC *m_pLArTPC;
    float m_tpcMinXEdge;
    float m_tpcMaxXEdge;
    unsigned int m_maxLoopIterations;
    float   m_distanceToLine;
    float m_boundaryTolerance;
    
};




//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

    
    
} // namespace lar_content

#endif // #ifndef LAR_TRACK_EXTENSION_REFINEMENT_ALGORITHM_H
