/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H
#define LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H 1

//#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRefinementBaseAlgorithm.h"

namespace lar_content
{
    
class CosmicRayEndpointCorrectionAlgorithm : public CosmicRayTrackRefinementBaseAlgorithm<ClusterEndpointAssociation>
{
public:
    
    CosmicRayEndpointCorrectionAlgorithm();
    
private:

    //pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterAssociationVector &clusterAssociationVector);
    
    bool IsDeltaRay(const pandora::Cluster *const pCluster, pandora::CartesianVector &clusterMergePoint, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream) const;

    void CreateMainTrack(ClusterEndpointAssociation &clusterEndpointAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterAssociationVector &clusterAssociationVector) const;

    
    void UpdateAfterMainTrackModification(const pandora::Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterAssociationVector &clusterAssociationVector) const;

    void RemoveClusterAssociationFromClusterVector(const ClusterEndpointAssociation &clusterAssociation, pandora::ClusterVector &clusterVector) const;

    void GetExtrapolatedCaloHits(ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, ClusterToCaloHitListMap &clusterToCaloHitListMap) const;

    bool IsExtrapolatedEndpointNearBoundary(const ClusterEndpointAssociation &clusterAssociation, const float boundaryTolerance) const;

    int m_minCaloHits;
    float m_maxDistanceFromTPC;
    float m_minScaledZOffset;
    float m_thresholdAngleDeviation;
    float m_thresholdAngleDeviationBetweenLayers;
    int m_maxAnomalousPoints;
    float m_thresholdMaxAngleDeviation;
    int m_deltaRayslidingFitWindow;
    float m_growingFitInitialLength;
    float m_growingFitSegmentLength;
    float m_furthestDistanceToLine;
    float m_closestDistanceToLine;
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H
