/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/SpanAcceptanceTool.h
 *
 *  @brief  Header file for the short span tool class
 *
 *  $Log: $
 */
#ifndef SPAN_ACCEPTANCE_TOOL_H
#define SPAN_ACCEPTANCE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  SpanAcceptanceTool class
 */
class SpanAcceptanceTool : public DeltaRayTensorTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    SpanAcceptanceTool();

    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void InvestigateSpanAcceptances(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade);
    bool IsConnected(const TensorType::Element &element) const;

};

} // namespace lar_content

#endif // #ifndef SHORT_SPAN_TOOL_H
