/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatching.h
 *
 *  @brief  Header file for the two view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{
    
class DeltaRayMatrixTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewDeltaRayMatchingAlgorithm class
 */
class TwoViewDeltaRayMatchingAlgorithm : public NViewMatchingAlgorithm<TwoViewMatchingControl<TrackTwoViewTopologyOverlapResult> >
{
public:
    typedef NViewMatchingAlgorithm<TwoViewMatchingControl<TrackTwoViewTopologyOverlapResult> > BaseAlgorithm;
    typedef std::map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterProximityMap;    
    typedef std::map<const pandora::Cluster*, pandora::CaloHitList> HitOwnershipMap;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef TwoViewDeltaRayMatchingAlgorithm::MatchingType::MatrixType MatrixType;
    typedef std::vector<pandora::HitType> HitTypeVector;
    
    /**
     *  @brief  Default constructor
     */
    TwoViewDeltaRayMatchingAlgorithm();

    void GetConnectedElements(const pandora::Cluster *const pClusterA, const bool hasAssociatedMuon, MatrixType::ElementList &elementList, pandora::ClusterSet &checkedClusters);

    void GetUnambiguousElements(const bool hasAssociatedMuon, MatrixType::ElementList &elementList);

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

    bool DoesClusterPassTesorThreshold(const pandora::Cluster *const pCluster) const;

    void RemoveThirdViewCluster(const pandora::Cluster *const pCluster);

    const ClusterToPfoMap &GetClusterToPfoMap(const pandora::HitType &hitType);

    const std::string &GetThirdViewClusterListName() const;

    pandora::StatusCode PerformThreeViewMatching(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
        float &reducedChiSquared) const;
    
    pandora::StatusCode PerformThreeViewMatching(const pandora::CaloHitList &clusterU, const pandora::CaloHitList &clusterV, const pandora::CaloHitList &clusterW,
        float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &XOverlap) const;

    bool CreatePfo(const MatrixType::Element &element);
    
private:
    void PrepareInputClusters(pandora::ClusterList &selectedClusters);

    void FillHitToClusterMap(const pandora::HitType &hitType);
    void AddToClusterMap(const pandora::Cluster *const pCluster);
    void FillClusterProximityMap(const pandora::HitType &hitType);
    void BuildKDTree(const pandora::HitType &hitType);
    void AddToClusterProximityMap(const pandora::Cluster *const pCluster);
    void FillClusterToPfoMap(const pandora::HitType &hitType);
    void FillStrayClusterList(const pandora::HitType &hitType);

    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3);
    
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        TrackTwoViewTopologyOverlapResult &overlapResult) const;

    void FindCommonMuonParents(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, pandora::PfoList &commonMuonPfoList) const;

    void GetNearbyMuonPfos(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::PfoList &nearbyMuonPfos) const;    

    pandora::StatusCode GetProjectedPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianPointVector &projectedPositions) const;
    
    void CollectThirdViewClusters(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::CartesianPointVector &projectedPositions,
        pandora::ClusterList &matchedClusters) const;

    void GetBestMatchedCluster(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::PfoList &commonMuonPfoList,
        const pandora::ClusterList &matchedClusters, const pandora::Cluster *&pBestMatchedCluster, float &reducedChiSquared) const;
    
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
        const pandora::CaloHitList &projectedHits, TrackTwoViewTopologyOverlapResult &overlapResult) const;

    void ExamineOverlapContainer();
    void TidyUp();    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MergeThirdView(const MatrixType::Element &element, const pandora::Cluster *const pSeedCluster);

    pandora::StatusCode PullOutDeltaRayHits(const MatrixType::Element &element,
        const pandora::ParticleFlowObject *const pParentMuon, pandora::CaloHitList &collectedHits) const;

    pandora::StatusCode ProjectMuonPositions(const pandora::HitType &thirdViewHitType, const pandora::ParticleFlowObject *const pParentMuon,
        pandora::CartesianPointVector &projectedPositions) const;

    pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &referencePoint, const pandora::CartesianPointVector &cartesianPointVector,
        const pandora::Cluster *const pCluster) const;

    float GetClosestDistance(const pandora::Cluster *const pCluster, const pandora::CartesianPointVector &cartesianPointVector) const;

    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &cartesianPointVector) const;

    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList) const;

    void GetClosestPositions(const pandora::CartesianPointVector &pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianVector &outputPosition1,
        pandora::CartesianVector &outputPosition2) const;

    const pandora::Cluster *SplitCluster(const pandora::Cluster *const pMuonCluster, pandora::CaloHitList &collectedHits) const;

    void GetBestMatchedAvailableCluster(const pandora::ClusterList &matchedClusters, const pandora::Cluster *&pBestMatchedCluster) const;

    void GrowThirdView(const MatrixType::Element &element, ProtoParticle &protoParticle);

    
    void GetClusterSpanX(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax) const;
    pandora::StatusCode GetClusterSpanZ(const pandora::CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const;

    //WILL NEED AN UPDATE AND DELETE CLUSTER FUNCTION

    HitToClusterMap       m_hitToClusterMap1;
    HitToClusterMap       m_hitToClusterMap2;

    HitKDTree2D m_kdTree1;    
    HitKDTree2D m_kdTree2;

    ClusterProximityMap   m_clusterProximityMap1;
    ClusterProximityMap   m_clusterProximityMap2;
    
    ClusterToPfoMap       m_clusterToPfoMap1;
    ClusterToPfoMap       m_clusterToPfoMap2;

    pandora::ClusterList m_strayClusterList1;
    pandora::ClusterList m_strayClusterList2;
    
    typedef std::vector<DeltaRayMatrixTool*> MatrixToolVector;
    MatrixToolVector                  m_algorithmToolVector;          ///< The algorithm tool vector
    
    unsigned int                      m_nMaxMatrixToolRepeats;
    unsigned int                      m_minClusterCaloHits;
    float                             m_searchRegion1D;             ///< Search region, applied to each dimension, for look-up from kd-tree        
    float                             m_xOverlapWindow;             ///< The maximum allowed displacement in x position
    float                             m_maxDisplacementSquared;
    float                             m_minMatchedFraction;
    unsigned int                      m_minMatchedPoints;
    unsigned int                      m_minProjectedPositions;
    float                      m_maxDistanceFromPrediction;
    float                      m_maxGoodMatchReducedChiSquared;
    float m_pseudoChi2Cut;
    std::string m_inputClusterListName;
    std::string m_muonPfoListName;
    std::string m_deltaRayPfoListName;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &TwoViewDeltaRayMatchingAlgorithm::GetThirdViewClusterListName() const
{
    return m_inputClusterListName;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoViewDeltaRayMatchingAlgorithm::ClusterToPfoMap &TwoViewDeltaRayMatchingAlgorithm::GetClusterToPfoMap(const pandora::HitType &hitType)
{
    const unsigned int hitTypeIndex(m_matchingControl.m_hitTypeToIndexMap.at(hitType));      

    return (hitTypeIndex == 1) ? m_clusterToPfoMap1 : m_clusterToPfoMap2;
}

    
/**
 *  @brief  DeltaRayTensorTool class
 */
class DeltaRayMatrixTool : public pandora::AlgorithmTool
{
public:
    typedef TwoViewDeltaRayMatchingAlgorithm::MatchingType::MatrixType MatrixType;
    typedef std::vector<MatrixType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &matrixTensor) = 0;
    };
    
} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
