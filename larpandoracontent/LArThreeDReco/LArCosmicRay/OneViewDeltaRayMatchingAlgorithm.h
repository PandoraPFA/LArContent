/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/OneViewDeltaRayMatchingAlgorithm.h
 *
 *  @brief  Header file for the one viw delta ray matching algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_ONE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_ONE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingContainers.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{

/**
 *  @brief  OneViewDeltaRayMatchingAlgorithm class
 */
class OneViewDeltaRayMatchingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    OneViewDeltaRayMatchingAlgorithm();
    
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the input cluster list of a given hit type
     *
     *  @param  hitType the hit type of list to retrieve
     *
     *  @return  the input cluster list of the specified hit type
     */
    const pandora::ClusterList GetInputClusterList(const pandora::HitType hitType);

    /**
     *  @brief  Get the input cosmic ray pfo list
     *
     *  @return  the input cosmic ray pfo list
     */
    const pandora::PfoList GetMuonPfoList();

    /**
     *  @brief  Get the input delta ray pfo list
     *
     *  @return  the input delta ray pfo list
     */
    const pandora::PfoList GetDeltaRayPfoList();

    /**
     *  @brief  Use nearby muon pfos to project into other views and attempt to match the remaining delta ray clusters 
     *
     *  @param  hitType the hit type of the input cluster and of the map to add to
     */
    void PerformOneViewMatching(const pandora::HitType hitType);

    /**
     *  @brief  Determine whether an input cluster belongs to a cosmic ray pfo
     *
     *  @param  pCluster the address of the input cluster
     *
     *  @return whether the input cluster belongs to a cosmic ray pfo
     */
    bool IsMuonPfo(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Use nearby muon pfos to project into other views and attempt to add a remaining delta ray cluster into an existing delta ray pfo
     *
     *  @param  pAvailableCluster the address of the delta ray cluster to add
     *  @param  nearbyMuonPfoVector the vector of nearby cosmic ray pfos
     *
     *  @return whether the input cluster was aded into an existing delta ray pfo
     */
    bool AddIntoExistingDeltaRay(const pandora::Cluster *const pAvailableCluster, const pandora::PfoVector &nearbyMuonPfoVector);

    /**
     *  @brief  Get the best matched available or unavailable cluster of a remaining delta ray cluster group wrt a cosmic ray pfo
     *
     *  @param  deltaRayClusterGroup the input group of to be merged remaining delta ray clusters
     *  @param  pNearbyMuonPfo the address of the nearby cosmic ray pfo
     *  @param  hitType the hit type to project into
     *  @param  findAvailable whether to search for available or unavailable delta ray clusters
     *
     *  @return whether the best matched cluster of the specified view
     */
    const pandora::Cluster *GetBestProjectedCluster(const pandora::ClusterList &deltaRayClusterGroup, const pandora::ParticleFlowObject *const pNearbyMuonPfo,
        const pandora::HitType hitType, const bool findAvailable);

    /**
     *  @brief  Determine whether an input cluster belongs to a delta ray pfo
     *
     *  @param  pCluster the address of the input cluster
     *
     *  @return whether the input cluster belongs to a delta ray pfo
     */
    bool IsDeltaRayPfo(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Determine cluster span (in x) of a group of clusters 
     *
     *  @param  clusterList the input list of clusters
     *  @param  spanMinX the output minimum x value
     *  @param  spanMaxX the output maximum x value
     */
    void GetClusterSpanX(const pandora::ClusterList &clusterList, float &spanMinX, float &spanMaxX);

    /**
     *  @brief  Use nearby muon pfos to project into other views and attempt to match a remaining delta ray cluster to form a delta ray pfo
     *
     *  @param  pAvailableCluster the address of the delta ray cluster to add
     *  @param  nearbyMuonPfoVector the vector of nearby cosmic ray pfos
     *  @param  modifiedClusters the output list of any delta ray clusters in the view of pAvailableCluster made unavailable in this process
     */
    void CreateDeltaRay(const pandora::Cluster *const pAvailableCluster, const pandora::PfoVector &nearbyMuonPfoVector, pandora::ClusterSet &modifiedClusters);

    /**
     *  @brief  In the view of the input available cluster, gather nearby available clusters
     *
     *  @param  pCluster the input available cluster
     *  @param  consideredClusters the list of investigated clusters 
     *  @param  foundClusters the output list of nearby available clusters (including the available cluster)
     */
    void GetNearbyAvailableClusters(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::ClusterList &foundClusters);
    
    /**
     *  @brief  Merge a collection of available clusters together updating hit containers accordingly
     *
     *  @param  clusterGroup the input group of clusters to merge
     *
     *  @param  the merged cluster
     */
    const pandora::Cluster *MergeClusterGroup(const pandora::ClusterList &clusterGroup);

    /**
     *  @brief  Create a pfo from the input clusters updating the cluster to pfo map accordingly
     *
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  pCluster3 the address of the third cluster
     */
    void CreatePfo(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3);

    /**
     *  @brief  Create a delta ray pfo from any remaining, significant clusters
     *
     *  @param  hitType the view to run in
     */
    void PerformRecovery(const pandora::HitType hitType);

    std::string m_muonPfoListName; ///< The list of reconstructed cosmic ray pfos
    std::string m_deltaRayPfoListName; ///< The list of reconstructed delta ray pfos
    std::string m_inputClusterListNameU; ///< The list of reconstructed U clusters
    std::string m_inputClusterListNameV; ///< The list of reconstructed V clusters
    std::string m_inputClusterListNameW; ///< The list of reconstructed W clusters
    std::string m_outputPfoListName; ///< The list to receive the created delta ray pfos
    DeltaRayMatchingContainers m_deltaRayMatchingContainers; ///< The class of hit, cluster and pfo ownership and proximity maps
    float m_overlapExtension; ///< The extension to each side of the x overlap region in which to search for matched clusters
    unsigned int m_minClusterHits; ///< The minimum number of hits for a cluster to be significant
};

} // namespace lar_content

#endif // #ifndef LAR_ONE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
