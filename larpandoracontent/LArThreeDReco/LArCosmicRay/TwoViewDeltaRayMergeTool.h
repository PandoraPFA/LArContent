/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMergeTool.h
 *
 *  @brief  Header file for the two view delta ray merge tool class
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_DELTA_RAY_MERGE_TOOL_H
#define TWO_VIEW_DELTA_RAY_MERGE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayMergeTool class
 */
class TwoViewDeltaRayMergeTool : public DeltaRayMatrixTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    TwoViewDeltaRayMergeTool();

private:
    bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakeMerges(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix, bool &mergesMade) const;
    
    bool PickOutGoodMatches(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList) const;

    bool CreatePfo(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element) const;

    void MergeThirdView(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const pandora::Cluster *const pSeedCluster) const;

    pandora::StatusCode CollectDeltaRayHits(const MatrixType::Element &element, const pandora::CartesianPointVector &deltaRayProjectedPositions,
        const pandora::ParticleFlowObject *const pParentMuon, pandora::CaloHitList &collectedHits) const;

    void ProjectMuonPositions(const pandora::HitType &thirdViewHitType, const pandora::ParticleFlowObject *const pParentMuon,
        pandora::CartesianPointVector &projectedPositions) const;

    pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &referencePoint, const pandora::CartesianPointVector &cartesianPointVector,
        const pandora::Cluster *const pCluster) const;

    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &cartesianPointVector) const;

    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList) const;

    void GetClosestPositions(const pandora::CartesianPointVector &pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianVector &outputPosition1,
        pandora::CartesianVector &outputPosition2) const;

    void ProjectPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianPointVector &projectedPositions) const;

    const pandora::Cluster *SplitCluster(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const pandora::Cluster *const pMuonCluster, pandora::CaloHitList &collectedHits) const;

    void GetBestMatchedAvailableCluster(const pandora::ClusterList &matchedClusters, const pandora::Cluster *&pBestMatchedCluster) const;

    void GrowThirdView(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, ProtoParticle &protoParticle) const;
    
    float m_maxUnambiguousClusterSeparation;
    float m_maxDRSeparationFromTrack;
    float m_maxVertexSeparation;
    float m_maxClusterSeparation;
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_DELTA_RAY_MERGE_TOOL_H
