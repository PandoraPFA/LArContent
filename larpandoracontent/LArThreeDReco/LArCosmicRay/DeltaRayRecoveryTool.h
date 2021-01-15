/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRecoveryTool.h
 *
 *  @brief  Header file for the delta ray recovery tool class
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_RECOVERY_TOOL_H
#define DELTA_RAY_RECOVERY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayRecoveryTool class
 */
class DeltaRayRecoveryTool : public DeltaRayTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayRecoveryTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakeMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor, bool &mergesMade) const;
    
    bool PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const;
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_RECOVERY_TOOL_H
