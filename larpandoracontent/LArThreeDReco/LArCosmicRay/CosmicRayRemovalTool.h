/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayRemoval.h
 *
 *  @brief  Header file for the long span tool class
 *
 *  $Log: $
 */
#ifndef COSMIC_RAY_REMOVAL_TOOL_H
#define COSMIC_RAY_REMOVAL_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayRemovalTool class
 */
class CosmicRayRemovalTool : public DeltaRayTensorTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    CosmicRayRemovalTool();

    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void SearchForMuonContamination(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade);

    bool PassElementChecks(TensorType::Element &element, const pandora::HitType &hitType);

    bool ShouldSplitDeltaRay(const pandora::Cluster *const pMuonCluster, const pandora::Cluster *const pDeltaRayCluster, const pandora::CartesianVector &muonPosition,
        const TwoDSlidingFitResult &slidingFitResult) const;

    void FindExtrapolatedHits(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary,
			      pandora::CaloHitList &collectedHits) const;

    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineEnd, const float distanceToLine) const;   

    void CreateSeed(const TensorType::Element &element, const pandora::HitType &badHitType, pandora::CaloHitList &collectedHits) const;

    bool GrowSeed(const TensorType::Element &element, const pandora::HitType &badHitType, pandora::CaloHitList &collectedHits) const;
    
    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList) const;

    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &cartesianPointVector) const;
    void ProjectPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianPointVector &projectedPositions) const;

    void SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &badHitType, pandora::CaloHitList &collectedHits) const;

    pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &referencePoint, const pandora::CartesianPointVector &cartesianPointVector, const pandora::Cluster *const pCluster) const;

    float m_minViewXSpan;
    float m_minSeparation;
    float m_xOverlapWindow;
};

} // namespace lar_content

#endif // #ifndef COSMIC_RAY_REMOVAL_TOOL_H
