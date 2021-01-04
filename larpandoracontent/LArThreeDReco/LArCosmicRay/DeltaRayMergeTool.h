/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMergeTool.h
 *
 *  @brief  Header file for the delta ray merge tool class
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_MERGE_TOOL_H
#define DELTA_RAY_MERGE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayMergeTool class
 */
class DeltaRayMergeTool : public DeltaRayTensorTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    DeltaRayMergeTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakeMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor, bool &mergesMade) const;
    
    bool MakeTwoCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, pandora::ClusterSet &modifiedClusters) const;

    void CombineCommonMuonPfoLists(const pandora::PfoList &commonMuonPfoList1, const pandora::PfoList &commonMuonPfoList2, pandora::PfoList &commonMuonPfoList) const;
    
    bool AreAssociated(const pandora::PfoList &commonMuonPfoList, const pandora::Cluster *const ClusterToEnlarge, const pandora::Cluster *const pClusterToDelete) const;

    bool IsHiddenTrack(const pandora::PfoList &commonMuonPfoList, const pandora::Cluster *const pClusterToEnlarge, const pandora::Cluster *const pClusterToDelete, bool &areAttached) const;

    bool IsConnected(const pandora::Pfo *const pCommonMuonPfo, const pandora::Cluster *const pCluster) const;

    void FindVertices(const pandora::Pfo *const pCommonMuonPfo, const pandora::Cluster *const pCluster, pandora::CaloHitList &vertexList) const;
    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList) const;

    bool IsBrokenCluster(const pandora::Cluster *const pClusterToEnlarge, const pandora::Cluster *const pClusterToDelete) const;

    bool MakeOneCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, pandora::ClusterSet &modifiedClusters) const;

    void PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, pandora::ClusterSet &modifiedClusters) const;

    float m_maxUnambiguousClusterSeparation;
    float m_maxDRSeparationFromTrack;
    float m_maxVertexSeparation;
    float m_maxClusterSeparation;
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_MERGE_TOOL_H
