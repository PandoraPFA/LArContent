/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatching.h
 *
 *  @brief  Header file for the three view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{
    
class DeltaRayTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewDeltaRayMatchingAlgorithm class
 */
class ThreeViewDeltaRayMatchingAlgorithm : public NViewMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult> >
{
public:
    typedef NViewMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult> > BaseAlgorithm;
    typedef std::map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterProximityMap;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef ThreeViewDeltaRayMatchingAlgorithm::MatchingType::TensorType TensorType;

    /**
     *  @brief  Default constructor
     */
    ThreeViewDeltaRayMatchingAlgorithm();

    void GetConnectedElements(const pandora::Cluster *const pCluster1, const bool hasAssociatedMuon, TensorType::ElementList &elementList, pandora::ClusterSet &checkedClusters);

    void GetUnambiguousElements(const bool hasAssociatedMuon, TensorType::ElementList &elementList);

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

    bool DoesClusterPassTesorThreshold(const pandora::Cluster *const pCluster) const;    

    void PrepareInputClusters(pandora::ClusterList &preparedClusterList);

    pandora::StatusCode PerformMatching(const pandora::CaloHitList &clusterU, const pandora::CaloHitList &clusterV, const pandora::CaloHitList &clusterW,
        float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &XOverlap) const;
    
    void UpdateForNewClusters(const pandora::ClusterVector &newClusterList, const pandora::PfoVector &pfoList);
    
    void UpdateContainers(const pandora::ClusterVector &newClusterVector, const pandora::PfoVector &pfoVector);
    
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);

    bool CreatePfos(ProtoParticleVector &protoParticleVector);

    std::string GetClusteringAlgName() const;
    
private:

    void FillHitToClusterMap(const pandora::HitType &hitType);
    void AddToClusterMap(const pandora::Cluster *const pCluster); 
    void FillClusterProximityMap(const pandora::HitType &hitType);
    void BuildKDTree(const pandora::HitType &hitType);
    void AddToClusterProximityMap(const pandora::Cluster *const pCluster);    
    void FillClusterToPfoMap(const pandora::HitType &hitType);
    void FillStrayClusterList(const pandora::HitType &hitType);
    
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        DeltaRayOverlapResult &overlapResult) const;

    void FindCommonMuonParents(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW, pandora::PfoList &commonMuonPfoList) const;

    void GetNearbyMuonPfos(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::PfoList &nearbyMuonPfos) const;

    void GetClusterSpanX(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax) const;
    pandora::StatusCode GetClusterSpanZ(const pandora::CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const;

    void CollectStrayHits(const pandora::Cluster *const pBadCluster, const float spanMinX, const float spanMaxX, pandora::ClusterList &collectedClusters);

    const pandora::ClusterList &GetStrayClusterList(const pandora::HitType &hitType) const;

    void AddInStrayClusters(const pandora::Cluster *const pClusterToEnlarge, const pandora::ClusterList &collectedClusters);

    void ExamineOverlapContainer();
    
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
    
    typedef std::vector<DeltaRayTensorTool*> TensorToolVector;
    TensorToolVector                  m_algorithmToolVector;          ///< The algorithm tool vector
    
    unsigned int                      m_nMaxTensorToolRepeats;
    unsigned int                      m_minClusterCaloHits;
    float                             m_searchRegion1D;             ///< Search region, applied to each dimension, for look-up from kd-tree    
    float                             m_pseudoChi2Cut;              ///< Pseudo chi2 cut for three view matching 
    float                             m_xOverlapWindow;             ///< The maximum allowed displacement in x position
    float                             m_minMatchedFraction;
    unsigned int                      m_minMatchedPoints;

    std::string  m_reclusteringAlgorithmName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string ThreeViewDeltaRayMatchingAlgorithm::GetClusteringAlgName() const
{
    return m_reclusteringAlgorithmName;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
/**
 *  @brief  DeltaRayTensorTool class
 */
class DeltaRayTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewDeltaRayMatchingAlgorithm::MatchingType::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};
    
} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
