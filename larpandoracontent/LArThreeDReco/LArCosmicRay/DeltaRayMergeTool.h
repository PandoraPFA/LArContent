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
    /**
     *  @brief  Default constructor
     */
    DeltaRayMergeTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor, bool &mergesMade) const;
    
    bool MakeTwoCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const;

    void CombineCommonMuonPfoLists(const pandora::PfoList &commonMuonPfoList1, const pandora::PfoList &commonMuonPfoList2, pandora::PfoList &commonMuonPfoList) const;

    bool AreAssociated(const TensorType::Element &element1, const TensorType::Element &element2, const pandora::HitType &mergeHitType) const;

    void GetConnectedMuons(const pandora::PfoList &commonMuonPfoList, const pandora::Cluster *const pClusterToEnlarge, pandora::PfoList &connectedMuonPfoList) const;

    bool IsHiddenTrack(const pandora::ParticleFlowObject *const pMuonPfo, const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;
    
    bool IsConnected(const pandora::Pfo *const pCommonMuonPfo, const pandora::Cluster *const pCluster) const;

    void FindVertices(const pandora::Pfo *const pCommonMuonPfo, const pandora::Cluster *const pCluster, pandora::CaloHitList &vertexList) const;

    bool IsBrokenCluster(const pandora::Cluster *const pClusterToEnlarge, const pandora::Cluster *const pClusterToDelete) const;

    bool MakeOneCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const;

    float m_maxUnambiguousClusterSeparation;
    float m_maxDRSeparationFromTrack;
    float m_maxVertexSeparation;
    float m_maxClusterSeparation;
    float m_maxGoodMatchReducedChiSquared;
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_MERGE_TOOL_H
