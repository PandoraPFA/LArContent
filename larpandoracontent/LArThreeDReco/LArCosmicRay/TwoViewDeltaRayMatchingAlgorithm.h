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
    typedef std::unordered_map<const pandora::CaloHit*, pandora::CaloHitList> HitAssociationMap;
    typedef std::map<const pandora::Cluster*, pandora::CaloHitList> HitOwnershipMap;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    
    /**
     *  @brief  Default constructor
     */
    TwoViewDeltaRayMatchingAlgorithm();

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

private:
    void PrepareInputClusters(pandora::ClusterList &selectedClusters);

    void FillHitToClusterMap(const pandora::HitType &hitType);
    void FillClusterProximityMap(const pandora::HitType &hitType);
    void FillClusterToPfoMap(const pandora::HitType &hitType);
    void FillHitAssociationMap();    
    
    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3);
    
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const HitAssociationMap &hitAssociationMap, TrackTwoViewTopologyOverlapResult &overlapResult) const;

    void GetNearbyMuonPfos(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::PfoList &nearbyMuonPfos) const;
    
    bool AreClustersCompatible(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    pandora::StatusCode GetProjectedPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianPointVector &projectedPositions) const;

    void CollectHits(const pandora::CartesianPointVector &projectedPositions, const HitAssociationMap &hitAssociationMap,
        const float xMin, const float xMax, pandora::CaloHitList &collectedHits, HitOwnershipMap &hitOwnershipMap) const;
    
    void CollectAssociatedHits(const pandora::CaloHit *const pSeedCaloHit, const pandora::CaloHit *const pCurrentCaloHit,
        const HitAssociationMap &hitAssociationMap, const float xMin, const float xMax, pandora::CaloHitList &associatedHitList) const;

    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
        const pandora::CaloHitList &projectedHits, TrackTwoViewTopologyOverlapResult &overlapResult) const;  

    void ExamineOverlapContainer();
    void TidyUp();    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //WILL NEED AN UPDATE AND DELETE CLUSTER FUNCTION

    
    HitToClusterMap       m_hitToClusterMap1;
    HitToClusterMap       m_hitToClusterMap2;

    HitKDTree2D m_kdTree1;    
    HitKDTree2D m_kdTree2;

    ClusterProximityMap   m_clusterProximityMap1;
    ClusterProximityMap   m_clusterProximityMap2;
    
    ClusterToPfoMap       m_clusterToPfoMap1;
    ClusterToPfoMap       m_clusterToPfoMap2;
    
    typedef std::vector<DeltaRayMatrixTool*> MatrixToolVector;
    MatrixToolVector                  m_algorithmToolVector;          ///< The algorithm tool vector
    
    unsigned int                      m_nMaxMatrixToolRepeats;
    unsigned int                      m_minClusterCaloHits;
    float                             m_searchRegion1D;             ///< Search region, applied to each dimension, for look-up from kd-tree        
    float                             m_xOverlapWindow;             ///< The maximum allowed displacement in x position
    float                             m_maxDisplacementSquared;
    float                             m_minMatchedFraction;
    unsigned int                      m_minMatchedPoints;
    float m_pseudoChi2Cut;
    std::string m_inputClusterListName;
    std::string m_muonPfoListName;
    std::string m_deltaRayPfoListName;
    HitAssociationMap m_hitAssociationMap;

};

//------------------------------------------------------------------------------------------------------------------------------------------
    
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
