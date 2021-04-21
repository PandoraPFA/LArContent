/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewAmbiguousDeltaRayTool.h
 *
 *  @brief  Header file for the two view delta ray merge tool class
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_AMBIGUOUS_DELTA_RAY_TOOL_H
#define TWO_VIEW_AMBIGUOUS_DELTA_RAY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewAmbiguousDeltaRayTool class
 */
class TwoViewAmbiguousDeltaRayTool : public DeltaRayMatrixTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    TwoViewAmbiguousDeltaRayTool();

private:
    bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void ExamineConnectedElements(MatrixType &overlapMatrix) const;
    
    bool PickOutGoodMatches(const MatrixType::ElementList &elementList) const;
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_AMBIGUOUS_DELTA_RAY_TOOL_H
