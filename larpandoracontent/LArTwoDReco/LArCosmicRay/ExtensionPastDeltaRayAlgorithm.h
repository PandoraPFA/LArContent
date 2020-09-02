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
 *  @brief ExtensionPastDeltaRayAlgorithm class
 */
class ExtensionPastDeltaRayAlgorithm :  public TrackExtensionRefinementAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */        
    ExtensionPastDeltaRayAlgorithm();
    
private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Whether the cluster endpoint meets the criteria: does the endpoint lie at the end of a clustering error
     *
     *  @param  microSlidingFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  clusterMergeDirection the merge direction of the cluster
     *  @param  isEndUpstream whether the cluster endpoint is at the upstream end of the cluster (i.e. lower z component)
     *  @param  pClusterList the list of all clusters
     *  @param  clusterMergePoint the merge position of the cluster
     *
     *  @return  whether a the cluster endpoint was found to lie at the end of a clustering error 
     */         
    bool DoesPassCriteria(const TwoDSlidingFitResult &microSlidingFitResult, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream, const pandora::ClusterList *const pClusterList, pandora::CartesianVector &clusterMergePoint) const;

    /**
     *  @brief  Whether a significant curve is found at the cluster merge point end of the hits or at the endpoint end of the hits
     *
     *  @param  microSlidingFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  clusterMergeDirection the merge direction of the cluster
     *  @param  isEndUpstream whether the cluster endpoint is at the upstream end of the cluster (i.e. lower z component)
     *  @param  clusterMergePoint the merge position of the cluster
     *
     *  @return  whether a curve was found
     */         
    bool IsDeltaRay(const TwoDSlidingFitResult &microSlidingFitResult, const pandora::CartesianVector &clusterMergeDirection, const bool isEndUpstream,
        pandora::CartesianVector &clusterMergePoint) const; 

    /**
     *  @brief  Whether a significant curve is found at a specified end of the hits between the cluster merge point and endpoint
     *
     *  @param  subsetFit the local TwoDSlidingFitResult of the hits in the ambiguous region
     *  @param  isEndUpstream whether the cluster endpoint is at the upstream end of the cluster (i.e. lower z component)
     *  @param  isClusterMergePointEnd whether the cluster merge point end of the ambiguous section is being investigated
     *  @param  clusterMergeDirection the merge direction of the cluster
     *  @param  clusterMergePoint the merge position of the cluster
     *
     *  @return  whether a curve was found
     */         
    bool IsCurvePresent(const TwoDSlidingFitResult &subsetFit, const bool isEndUpstream, const bool isClusterMergePointEnd,
        const pandora::CartesianVector &clusterMergeDirection, pandora::CartesianVector &clusterMergePoint) const;
    
    float           m_maxDistanceFromTPC;                     ///< The distance from the TPC boundary in which the clusterMergePoint or endpoint must fall
    float           m_minZOffset;                             ///< The threshold z deviation between the cluster endpoint and extrapolated endpoint of a clustering error
    unsigned int    m_deltaRayslidingFitWindow;               ///< The sliding fit window used to compute the TwoDSlidingFitResult of the hits in the ambiguous region
    float           m_thresholdAngleDeviation;                ///< The opening angle between the subsetFit layer and clusterMergeDirection that suggests the start of a curve
    float           m_thresholdMaxAngleDeviation;             ///< The required opening angle between the subsetFit layer and clusterMergeDirection that suggests the curve is significant
    float           m_thresholdAngleDeviationBetweenLayers;   ///< The amount the opening angle is required to increase between layers for the curve to be continuous
    int             m_maxAnomalousPoints;                     ///< The number of 'wobbles' the curve is allowed to make
};

} // namespace lar_content

#endif // #ifndef LAR_EXTENSION_PAST_DELTA_RAY_ALGORITHM_H
