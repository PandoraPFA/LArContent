/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/AmbiguousDeltaRayTool.h
 *
 *  @brief  Header file for the ambiguous delta ray tool class.
 *
 *  $Log: $
 */
#ifndef AMBIGUOUS_DELTA_RAY_TOOL_H
#define AMBIGUOUS_DELTA_RAY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  AmbiguousDeltaRayTool class
 */
class AmbiguousDeltaRayTool : public DeltaRayTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    AmbiguousDeltaRayTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) const;

    bool PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const;

    float m_maxGoodMatchReducedChiSquared;
};

} // namespace lar_content

#endif // #ifndef AMBIGUOUS_DELTA_RAY_TOOL_H
