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

    /**
     *  @brief  Identify ambiguous matches (e.g. 3:2:1) and, if possible, create pfos out of the best 1:1:1 cluster match
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     */
    void ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) const;

    /**
     *  @brief  Identify the best 1:1:1 match in a group of connected elements and from it create a pfo
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  elementList the tensor element list
     *
     *  @return whether any particles were created
     */
    bool PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const;

    float m_maxGoodMatchReducedChiSquared; ///< The maximum reduced chi squared value of a good 1:1:1 match
};

} // namespace lar_content

#endif // #ifndef AMBIGUOUS_DELTA_RAY_TOOL_H
