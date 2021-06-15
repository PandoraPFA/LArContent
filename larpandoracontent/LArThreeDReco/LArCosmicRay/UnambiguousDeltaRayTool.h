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

    /**
     *  @brief  Create delta ray pfos out of unambiguous (1:1:1) matches that are connected to a parent cosmic ray
     *
     *  @param  elementList the tensor element list
     *
     *  @return  whether any delta ray pfos were created
     */
    bool ExamineUnambiguousElements(TensorType::ElementList &elementList);

    /**
     *  @brief  Determine whether the clusters of an element are connected to a cosmic ray pfo
     *
     *  @param  elementList the tensor element
     *
     *  @return  whether the clusters are connected to a cosmic ray pfo
     */
    bool IsConnected(const TensorType::Element &element) const;

    float m_maxSeparation;                ///< The maximum separation between a connected delta ray cluster and a cosmic ray cluster
    unsigned int m_minNConnectedClusters; ///< The threshold number of connected delta ray clusters required for particle creation
};

} // namespace lar_content

#endif // #ifndef UNAMBIGUOUS_DELTA_RAY_TOOL
