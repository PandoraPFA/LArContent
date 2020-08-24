/*
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionThroughShowerAlgorithm.h
 *
 *  @brief  Header file for the extension through shower class.
 *
 *  $Log: $
 */
#ifndef LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H
#define LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/TrackExtensionRefinementAlgorithm.h"

namespace lar_content
{
/**
 *  @brief ExtensionThroughShowerAlgorithm class
 */
class ExtensionThroughShowerAlgorithm :  public TrackExtensionRefinementAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */    
    ExtensionThroughShowerAlgorithm();
    
private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  If it exists, find the cluster endpoint association of the furthest cluster from the tpc that is separated from the given TPC boundary by a shower 
     *
     *  @param  clusterVector the vector of clusters to consider
     *  @param  slidingFitResultMapPair the {micro, macro} pair of [cluster -> TwoDSlidingFitResult] maps 
     *  @param  pClusterList the list of all clusters
     *  @param  isHigherXboundary whether to look for endpoints that are closer to the higher or lower x tpc boundary
     *  @param  clusterAssociation the output cluster endpoint association
     *
     *  @return  bool whether a cluster endpoint association was found
     */        
    bool FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        const pandora::ClusterList *const pClusterList, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation);

    /**
     *  @brief  Whether the cluster endpoint lies in front a shower region 
     *
     *  @param  microSlidingFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  clusterMergePoint the merge position of the cluster
     *  @param  clusterMergeDirection the merge direction of the cluster
     *  @param  isEndUpstream whether the cluster endpoint is at the upstream end of the cluster (i.e. lower z component)
     *  @param  pClusterList the list of all clusters
     *
     *  @return  bool whether a the cluster endpoint was found to lie in front of a shower
     */         
    bool IsClusterEndpointInFrontOfShower(const TwoDSlidingFitResult &microSlidingFitResult, const pandora::CartesianVector &clusterMergePoint,
        const pandora::CartesianVector &clusterMergeDirection, const bool isEndUpstream, const pandora::ClusterList *const pClusterList);   

    float          m_maxTrackDistanceToShowerBranch;         ///< The maximum distance of a track cluster from a shower branch
    float          m_maxShowerBranchTransverseDistance;      ///< The maximum transverse projection of the closest pseudolayer position of a shower branch cluster  
    float          m_maxShowerBranchLongitudinalDistance;    ///< The maximum longitudinal projection of the closest pseudolayer position of a shower branch cluster  
    unsigned int   m_thresholdShowerClusterCount;            ///< The threshold number of shower branches of a significant shower region
};
 
} // namespace lar_content

#endif // #ifndef LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H
