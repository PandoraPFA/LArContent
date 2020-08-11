/*
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionPastDeltaRayAlgorithm.h
 *
 *  @brief  Header file for the extension past delta ray class.
 *
 *  $Log: $
 */
#ifndef LAR_EXTENSION_PAST_DELTA_RAY_ALGORITHM_H
#define LAR_EXTENSION_PAST_DELTA_RAY_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackExtensionRefinementAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  class
 */
class ExtensionPastDeltaRayAlgorithm :  public TrackExtensionRefinementAlgorithm
{
public:
    
    ExtensionPastDeltaRayAlgorithm();
    
protected:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterEndpointAssociation &clusterAssociation, const pandora::ClusterList *const pClusterList, const bool isHigherXBoundary);

    
    bool IsDeltaRay(const pandora::Cluster *const pCluster, pandora::CartesianVector &clusterMergePoint, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream) const;

    float m_maxDistanceFromTPC;
    float m_minScaledZOffset;
    float m_thresholdAngleDeviation;
    float m_thresholdAngleDeviationBetweenLayers;
    int m_maxAnomalousPoints;
    float m_thresholdMaxAngleDeviation;
    int m_deltaRayslidingFitWindow;
};
 
} // namespace lar_content

#endif // #ifndef LAR_EXTENSION_PAST_DELTA_RAY_ALGORITHM_H
