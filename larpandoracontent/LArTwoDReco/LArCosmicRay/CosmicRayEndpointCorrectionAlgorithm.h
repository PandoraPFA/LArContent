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

class CosmicRayEndpointCorrectionAlgorithm : public CosmicRayTrackRefinementBaseAlgorithm
{
public:
    
    CosmicRayEndpointCorrectionAlgorithm();

private:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    /**
     *  @brief  Select clusters to be considered in algorithm
     *
     *  @param  pClusterList the input cluster list
     *  @param  clusterVector the output cluster vector
     */
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector) const;

    bool IsDeltaRay(const ClusterAssociation &clusterAssociation, const bool isUpstream) const;

    
    bool InvestigateClusterEnd(const pandora::Cluster *const pCluster, const bool isUpstream, const TwoDSlidingFitResult &microFitResult, const TwoDSlidingFitResult &macroFitResult, ClusterAssociation &clusterAssociation);
    

    void RefineTrackEndpoint(const pandora::Cluster *const pCluster, const pandora::CartesianVector &clusterUpstreamMergePoint, const pandora::CartesianVector &clusterDownstreamMergePoint, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapVector) const;

    
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
