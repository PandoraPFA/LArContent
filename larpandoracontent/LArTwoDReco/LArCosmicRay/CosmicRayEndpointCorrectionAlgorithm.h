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
    
    bool IsDeltaRay(const pandora::Cluster *const pCluster, const pandora::CartesianVector &clusterMergePoint, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream) const;   

    void CreateMainTrack(ClusterEndpointAssociation &clusterEndpointAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    
    void UpdateAfterMainTrackModification(const pandora::Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPai) const;   

    
    int m_minCaloHits;
    float m_maxDistanceFromTPC;
    float m_curveThreshold;
    float m_minScaledZOffset;
    float m_thresholdAngleDeviation;
    float m_thresholdAngleDeviationBetweenLayers;
    int m_maxAnomalousPoints;
    float m_thresholdMaxAngleDeviation;
    
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H
