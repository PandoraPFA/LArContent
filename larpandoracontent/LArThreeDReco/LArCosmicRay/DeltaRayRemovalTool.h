/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.h
 *
 *  @brief  Header file for the delta ray removal tool class.
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_REMOVAL_TOOL_H
#define DELTA_RAY_REMOVAL_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"


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

    void ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
        bool &changesMade) const;

    bool PassElementChecks(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &hitType) const;

    bool IsContaminated(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &hitType) const;

    void SplitMuonCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
        const pandora::HitType &hitType, const pandora::CaloHitList &deltaRayHits) const;

    float m_minSeparation;
    unsigned int m_slidingFitWindow;
    float m_minDeviationFromTransverse;
    unsigned int m_significantHitThreshold;
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_REMOVAL_TOOL_H
