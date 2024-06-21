/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingContainers.h
 *
 *  @brief  Header file for the delta ray matching containers class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_MATCHING_CONTAINERS_H
#define LAR_DELTA_RAY_MATCHING_CONTAINERS_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayMatchingContainers class
 */
class DeltaRayMatchingContainers
{
public:
    typedef std::map<const pandora::Cluster *, const pandora::ParticleFlowObject *> ClusterToPfoMap;
    typedef std::map<const pandora::Cluster *, pandora::ClusterList> ClusterProximityMap;

    /**
     *  @brief  Default constructor
     */
    DeltaRayMatchingContainers();

    /**
     *  @brief  Get the mapping of clusters to the pfos to which they belong
     *
     *  @return  the cluster to pfo map
     */
    const ClusterToPfoMap &GetClusterToPfoMap(const pandora::HitType hitType) const;

    /**
     *  @brief  Get the mapping of clusters to to their neighbouring clusters
     *
     *  @return  the cluster to nearby cluster map
     */
    const ClusterProximityMap &GetClusterProximityMap(const pandora::HitType hitType) const;

    /**
     *  @brief  Fill the HitToClusterMap, the ClusterProximityMap and the ClusterToPfoMap in all input views
     *
     *  @param  inputPfoList the input list of pfos
     *  @param  inputClusterList1 the input list of clusters in view 1
     *  @param  inputClusterList2 the input list of clusters in view 2
     *  @param  inputClusterList3 the input list of clusters in view 3
     */
    void FillContainers(const pandora::PfoList &inputPfoList, const pandora::ClusterList &inputClusterList1,
        const pandora::ClusterList &inputClusterList2 = pandora::ClusterList(),
        const pandora::ClusterList &inputClusterList3 = pandora::ClusterList());

    /**
     *  @brief  Add the clusters of a cosmic ray/delta ray pfo to the cluster to pfo maps
     *
     *  @param  the address of the input cosmic ray/delta ray pfo
     */
    void AddClustersToPfoMaps(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Add a list of clusters to the hit to cluster and cluster proximity maps and, if appropriate, to the cluster to pfo map
     *
     *  @param  newClusterVector the ordered cluster vector
     *  @param  pfoVector the matching ordered vector of pfos to which the clusters belong (nullptr if not applicable)
     */
    void AddClustersToContainers(const pandora::ClusterVector &newClusterVector, const pandora::PfoVector &pfoVector);

    /**
     *  @brief  Remove an input cluster's hits from the hit to cluster and cluster proximity maps and, if appropriate, from the cluster to pfo map
     *
     *  @param  pDeletedCluster the input cluster
     */
    void RemoveClusterFromContainers(const pandora::Cluster *const pDeletedCluster);

    /**
     *  @brief  Empty all algorithm containers
     */
    void ClearContainers();

    float m_searchRegion1D; ///< Search region, applied to each dimension, for look-up from kd-tree

private:
    typedef std::map<const pandora::CaloHit *, const pandora::Cluster *> HitToClusterMap;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    /**
     *  @brief  Populate the hit to cluster map from a list of clusters
     *
     *  @param  inputClusterList the input list of clusters
     */
    void FillHitToClusterMap(const pandora::ClusterList &inputClusterList);

    /**
     *  @brief  Add the hits of a given cluster to the hit to cluster map
     *
     *  @param  pCluster the address of the input cluster
     */
    void AddToClusterMap(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Populate all cluster to pfo maps from a list of particle flow objects
     *
     *  @param  pfoList the input list of pfos
     */
    void FillClusterToPfoMaps(const pandora::PfoList &pfoList);

    /**
     *  @brief  Populate the cluster proximity map from a list of clusters
     *
     *  @param  inputClusterList the input list of clusters
     */
    void FillClusterProximityMap(const pandora::ClusterList &inputClusterList);

    /**
     *  @brief  Build the KD tree
     *
     *  @param  hitType the hit type of the KD tree to build
     */
    void BuildKDTree(const pandora::HitType hitType);

    /**
     *  @brief  Add a cluster to the cluster proximity map
     *
     *  @param  pCluster the address of the input cluster
     */
    void AddToClusterProximityMap(const pandora::Cluster *const pCluster);

    HitToClusterMap m_hitToClusterMapU;         ///< The mapping of hits to the clusters to which they belong (in the U view)
    HitToClusterMap m_hitToClusterMapV;         ///< The mapping of hits to the clusters to which they belong (in the V view)
    HitToClusterMap m_hitToClusterMapW;         ///< The mapping of hits to the clusters to which they belong (in the W view)
    HitKDTree2D m_kdTreeU;                      ///< The KD tree (in the U view)
    HitKDTree2D m_kdTreeV;                      ///< The KD tree (in the V view)
    HitKDTree2D m_kdTreeW;                      ///< The KD tree (in the W view)
    ClusterProximityMap m_clusterProximityMapU; ///< The mapping of clusters to their neighbouring clusters (in the U view)
    ClusterProximityMap m_clusterProximityMapV; ///< The mapping of clusters to their neighbouring clusters (in the V view)
    ClusterProximityMap m_clusterProximityMapW; ///< The mapping of clusters to their neighbouring clusters (in the W view)
    ClusterToPfoMap m_clusterToPfoMapU;         ///< The mapping of cosmic ray U clusters to the cosmic ray pfos to which they belong
    ClusterToPfoMap m_clusterToPfoMapV;         ///< The mapping of cosmic ray V clusters to the cosmic ray pfos to which they belong
    ClusterToPfoMap m_clusterToPfoMapW;         ///< The mapping of cosmic ray W clusters to the cosmic ray pfos to which they belong
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const DeltaRayMatchingContainers::ClusterProximityMap &DeltaRayMatchingContainers::GetClusterProximityMap(const pandora::HitType hitType) const
{
    return ((hitType == pandora::TPC_VIEW_U)   ? m_clusterProximityMapU
            : (hitType == pandora::TPC_VIEW_V) ? m_clusterProximityMapV
                                               : m_clusterProximityMapW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const DeltaRayMatchingContainers::ClusterToPfoMap &DeltaRayMatchingContainers::GetClusterToPfoMap(const pandora::HitType hitType) const
{
    return ((hitType == pandora::TPC_VIEW_U)   ? m_clusterToPfoMapU
            : (hitType == pandora::TPC_VIEW_V) ? m_clusterToPfoMapV
                                               : m_clusterToPfoMapW);
}

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_MATCHING_CONTAINERS_H
