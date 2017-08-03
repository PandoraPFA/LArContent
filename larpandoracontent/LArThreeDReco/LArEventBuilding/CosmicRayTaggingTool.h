/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/CosmicRayTaggingTool.h
 *
 *  @brief  Header file for the cosmic-ray tagging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TAGGING_TOOL_H
#define LAR_COSMIC_RAY_TAGGING_TOOL_H 1

#include "larpandoracontent/LArUtility/ParentAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayTaggingTool class
 */
class CosmicRayTaggingTool : public CosmicRayTaggingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRayTaggingTool();

    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TAGGING_TOOL_H
