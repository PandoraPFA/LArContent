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
    typedef std::map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterProximityMap;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const pandora::ClusterList GetInputClusterList(const pandora::HitType &hitType);
    const pandora::PfoList GetMuonPfoList();
    void FillHitToClusterMap(const pandora::HitType &hitType);
    void AddToClusterMap(const pandora::Cluster *const pCluster);
    void FillClusterToPfoMap(const pandora::HitType &hitType);
    void FillClusterProximityMap(const pandora::HitType &hitType);
    void BuildKDTree(const pandora::HitType &hitType);
    void AddToClusterProximityMap(const pandora::Cluster *const pCluster);

    void PerformOneViewMatching(const pandora::HitType &hitType);

    void GetNearbyClusters(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::ClusterList &foundClusters);
    
    void GetProjectedNearbyClusters(const pandora::ClusterList &deltaRayClusterGroup, const pandora::ParticleFlowObject *const pMuonPfo,
        const pandora::HitType &hitType, pandora::ClusterList &foundClusters);

    void GetClusterSpanX(const pandora::ClusterList &clusterList, float &xMin, float &xMax);

    void CreatePfos(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3);

    void RemoveClusterFromProximityMaps(const pandora::Cluster *const pCluster);

    const pandora::Cluster *MergeClusterGroup(const pandora::ClusterList &clusterGroup);

    void PerformRecovery(const pandora::HitType &hitType);
    
    HitToClusterMap       m_hitToClusterMapU;
    HitToClusterMap       m_hitToClusterMapV;
    HitToClusterMap       m_hitToClusterMapW;

    ClusterToPfoMap       m_clusterToPfoMapU;
    ClusterToPfoMap       m_clusterToPfoMapV;
    ClusterToPfoMap       m_clusterToPfoMapW;
    
    HitKDTree2D m_kdTreeU;    
    HitKDTree2D m_kdTreeV;
    HitKDTree2D m_kdTreeW;

    ClusterProximityMap   m_muonProximityMapU;
    ClusterProximityMap   m_muonProximityMapV;
    ClusterProximityMap   m_muonProximityMapW;
   
    ClusterProximityMap   m_availableProximityMapU;
    ClusterProximityMap   m_availableProximityMapV;
    ClusterProximityMap   m_availableProximityMapW;

    float    m_availableSearchRegion1D;
    
    std::string m_muonPfoListName;    
    std::string m_inputClusterListNameU;
    std::string m_inputClusterListNameV;
    std::string m_inputClusterListNameW;
    std::string m_outputPfoListName; 
};

} // namespace lar_content

#endif // #ifndef LAR_ONE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
