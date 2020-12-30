/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemoval.h
 *
 *  @brief  Header file for the delta ray removal tool
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_REMOVAL_TOOL_H
#define DELTA_RAY_REMOVAL_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h"

namespace lar_content
{
/**
 *  @brief  DeltaRayRemovalTool class
 */
class DeltaRayRemovalTool : public RemovalBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayRemovalTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);    

    void RemoveDeltaRayHits(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
        bool &changesMade) const;

    bool PassElementChecks(const TensorType::Element &element, const pandora::HitType &hitType) const;

    bool IsMuonEndpoint(const TensorType::Element &element) const;    

    bool IsContaminated(const TensorType::Element &element, const pandora::HitType &hitType) const;

    pandora::StatusCode CollectDeltaRayHits(const TensorType::Element &element, const pandora::HitType &badHitType,
        pandora::CaloHitList &collectedHits) const;

    void SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &badHitType,
        pandora::CaloHitList &collectedHits) const;

    float m_minSeparation;
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_REMOVAL_TOOL_H
