/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayRecoveryTool.h
 *
 *  @brief  Header file for the two view delta ray merge tool class
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_DELTA_RAY_RECOVERY_TOOL_H
#define TWO_VIEW_DELTA_RAY_RECOVERY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayRecoveryTool class
 */
class TwoViewDeltaRayRecoveryTool : public DeltaRayMatrixTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewDeltaRayRecoveryTool();

private:
    bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakeMerges(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix, bool &mergesMade) const;
    
    bool PickOutGoodMatches(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList) const;

};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_DELTA_RAY_RECOVERY_TOOL_H
