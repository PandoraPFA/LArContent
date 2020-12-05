/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ShortSpan.h
 *
 *  @brief  Header file for the short span tool class
 *
 *  $Log: $
 */
#ifndef SHORT_SPAN_TOOL_H
#define SHORT_SPAN_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ShortSpanTool class
 */
class ShortSpanTool : public DeltaRayTensorTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    ShortSpanTool();

    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    void InvestigateShortSpans(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade);

    bool GetShortCluster(const TensorType::Element &element, const pandora::Cluster *&pShortCluster) const;

    void GetGoodXOverlapExtrema(const TensorType::Element &element, const pandora::HitType &badHitType, float &minX, float &maxX) const;

    float m_minXOverlapFraction;
};

} // namespace lar_content

#endif // #ifndef SHORT_SPAN_TOOL_H
