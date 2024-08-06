/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.h
 *
 *  @brief  Header file for the cheating cosmic-ray tagging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_COSMIC_RAY_TAGGING_TOOL_H
#define LAR_CHEATING_COSMIC_RAY_TAGGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingBaseTool.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingCosmicRayTaggingTool class
 */
class CheatingCosmicRayTaggingTool : public CosmicRayTaggingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingCosmicRayTaggingTool();

    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxCosmicRayFraction; ///< The maximum cosmic ray fraction for a pfo to be declared an ambiguous cosmic ray
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_COSMIC_RAY_TAGGING_TOOL_H
