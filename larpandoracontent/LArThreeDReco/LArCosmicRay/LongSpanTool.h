/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/LongSpan.h
 *
 *  @brief  Header file for the long span tool class
 *
 *  $Log: $
 */
#ifndef LONG_SPAN_TOOL_H
#define LONG_SPAN_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  LongSpanTool class
 */
class LongSpanTool : public DeltaRayTensorTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    LongSpanTool();

    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void InvestigateLongSpans(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade);

    bool GetLongCluster(const TensorType::Element &element, const pandora::Cluster *&pLongCluster) const;

    bool IsConnected(const TensorType::Element &element) const;
    
    void GetGoodXOverlapExtrema(const TensorType::Element &element, const pandora::HitType &badHitType, float &minX, float &maxX) const;

    bool IsMuonEndpoint(const TensorType::Element &element, const pandora::HitType &badHitType) const;

    bool ShouldSplitDeltaRay(const pandora::Cluster *const pMuonCluster, const pandora::Cluster *const pDeltaRayCluster, const pandora::CartesianVector &muonPosition,
        const TwoDSlidingFitResult &slidingFitResult) const;

    void FindExtrapolatedHits(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary,
			      pandora::CaloHitList &collectedHits) const;

    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineEnd, const float distanceToLine) const;   

    void CreateSeed(const TensorType::Element &element, const pandora::HitType &badHitType, pandora::CaloHitList &collectedHits) const;

    void GrowSeed(const TensorType::Element &element, const pandora::HitType &badHitType, pandora::CaloHitList &collectedHits) const;
    
    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList) const;

    float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &cartesianPointVector) const;
    void ProjectPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianPointVector &projectedPositions) const;

    void SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &badHitType, pandora::CaloHitList &collectedHits) const;

    float m_xOverlapWindow;
    float m_minXOverlapFraction;
};

} // namespace lar_content

#endif // #ifndef LONG_SPAN_TOOL_H
