/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/UnambiguousDeltaRayTool.h
 *
 *  @brief  Header file for the unambiguous delta ray tool class.
 *
 *  $Log: $
 */
#ifndef UNAMBIGUOUS_DELTA_RAY_TOOL_H
#define UNAMBIGUOUS_DELTA_RAY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  UnambiguousDeltaRayTool class
 */
class UnambiguousDeltaRayTool : public DeltaRayTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    UnambiguousDeltaRayTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    void ExamineUnambiguousElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade);
    
    bool IsConnected(const TensorType::Element &element) const;

    float m_minSeparation;
    unsigned int m_minNConnectedClusters;
};

} // namespace lar_content

#endif // #ifndef SHORT_SPAN_TOOL_H
