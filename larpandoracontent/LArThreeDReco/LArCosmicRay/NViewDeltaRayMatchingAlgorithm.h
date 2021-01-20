/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatching.h
 *
 *  @brief  Header file for the n view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_N_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_N_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{
/**
 *  @brief  NViewDeltaRayMatchingAlgorithm class
 */
template<typename T>    
class NViewDeltaRayMatchingAlgorithm : public NViewMatchingAlgorithm<T>
{
public:
    //typedef NViewMatchingAlgorithm<T> BaseAlgorithm;
    
    typedef std::map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterProximityMap;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::vector<pandora::HitType> HitTypeVector;

    /**
     *  @brief  Default constructor
     */
    NViewDeltaRayMatchingAlgorithm();
    
    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

    virtual bool DoesClusterPassTesorThreshold(const pandora::Cluster *const pCluster) const = 0;

    void PrepareInputClusters(pandora::ClusterList &preparedClusterList);

    pandora::StatusCode PerformThreeViewMatching(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
        float &reducedChiSquared) const;
    
    pandora::StatusCode PerformThreeViewMatching(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &XOverlap) const;

    pandora::StatusCode PerformThreeViewMatching(const pandora::CaloHitList &clusterU, const pandora::CaloHitList &clusterV, const pandora::CaloHitList &clusterW,
        float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &XOverlap) const;
    
    void UpdateForNewClusters(const pandora::ClusterVector &newClusterList, const pandora::PfoVector &pfoList);
    
    void UpdateContainers(const pandora::ClusterVector &newClusterVector, const pandora::PfoVector &pfoVector);
    
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);

    bool CreatePfos(ProtoParticleVector &protoParticleVector);
    
protected:
    void FillHitToClusterMap(const pandora::HitType &hitType);
    void AddToClusterMap(const pandora::Cluster *const pCluster); 
    void FillClusterProximityMap(const pandora::HitType &hitType);
    void BuildKDTree(const pandora::HitType &hitType);
    void AddToClusterProximityMap(const pandora::Cluster *const pCluster);    
    void FillClusterToPfoMap(const pandora::HitType &hitType);
    void FillStrayClusterList(const pandora::HitType &hitType);
    
    void GetNearbyMuonPfos(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::PfoList &nearbyMuonPfos) const;

    void GetClusterSpanX(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax) const;
    
    pandora::StatusCode GetClusterSpanZ(const pandora::CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const;

    void CollectStrayHits(const pandora::Cluster *const pBadCluster, const float spanMinX, const float spanMaxX, pandora::ClusterList &collectedClusters);

    const pandora::ClusterList &GetStrayClusterList(const pandora::HitType &hitType) const;

    void AddInStrayClusters(const pandora::Cluster *const pClusterToEnlarge, const pandora::ClusterList &collectedClusters);

    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string   m_muonPfoListName;
    
    HitToClusterMap       m_hitToClusterMapU;
    HitToClusterMap       m_hitToClusterMapV;
    HitToClusterMap       m_hitToClusterMapW;

    HitKDTree2D m_kdTreeU;    
    HitKDTree2D m_kdTreeV;
    HitKDTree2D m_kdTreeW;

    ClusterProximityMap   m_clusterProximityMapU;
    ClusterProximityMap   m_clusterProximityMapV;
    ClusterProximityMap   m_clusterProximityMapW;
    
    ClusterToPfoMap       m_clusterToPfoMapU;
    ClusterToPfoMap       m_clusterToPfoMapV;
    ClusterToPfoMap       m_clusterToPfoMapW;

    pandora::ClusterList m_strayClusterListU;
    pandora::ClusterList m_strayClusterListV;
    pandora::ClusterList m_strayClusterListW;    
    
    float                             m_searchRegion1D;             ///< Search region, applied to each dimension, for look-up from kd-tree    
    float                             m_pseudoChi2Cut;              ///< Pseudo chi2 cut for three view matching 
    float                             m_xOverlapWindow;             ///< The maximum allowed displacement in x position
    float                             m_minMatchedFraction;
    unsigned int                      m_minMatchedPoints;
};
    
} // namespace lar_content

#endif // #ifndef LAR_N_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
