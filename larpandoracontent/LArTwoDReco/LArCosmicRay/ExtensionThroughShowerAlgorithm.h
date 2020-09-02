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
     *  @brief  Whether the cluster endpoint meets the criteria: does the endpoint lies in front a shower region 
     *
     *  @param  microSlidingFitResult the local TwoDSlidingFitResult of the cluster
     *  @param  clusterMergeDirection the merge direction of the cluster
     *  @param  isEndUpstream whether the cluster endpoint is at the upstream end of the cluster (i.e. lower z component)
     *  @param  pClusterList the list of all clusters
     *  @param  clusterMergePoint the merge position of the cluster
     *
     *  @return  whether the cluster endpoint was found to lie in front of a shower
     */         
    bool DoesPassCriteria(const TwoDSlidingFitResult &microSlidingFitResult, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream, const pandora::ClusterList *const pClusterList, pandora::CartesianVector &clusterMergePoint) const;  
    
    float          m_maxTrackDistanceToShowerBranch;         ///< The maximum distance of a track cluster from a shower branch
    float          m_maxShowerBranchTransverseDistance;      ///< The maximum transverse projection of the closest pseudolayer position of a shower branch cluster  
    float          m_maxShowerBranchLongitudinalDistance;    ///< The maximum longitudinal projection of the closest pseudolayer position of a shower branch cluster  
    unsigned int   m_thresholdShowerClusterCount;            ///< The threshold number of shower branches of a significant shower region
};
 
} // namespace lar_content

#endif // #ifndef LAR_EXTENSION_THROUGH_SHOWER_ALGORITHM_H
