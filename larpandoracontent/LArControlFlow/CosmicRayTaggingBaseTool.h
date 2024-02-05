/**
 *  @file   larpandoracontent/LArControlFlow/CosmicRayTaggingBaseTool.h
 *
 *  @brief  Header file for the cosmic ray tagging tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TAGGING_BASE_TOOL_H
#define LAR_COSMIC_RAY_TAGGING_BASE_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class MasterAlgorithm;

/**
 *  @brief  CosmicRayTaggingBaseTool class
 */
class CosmicRayTaggingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Find the list of ambiguous pfos (could represent cosmic-ray muons or neutrinos)
     *
     *  @param  parentCosmicRayPfos the list of parent cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous pfos
     *  @param  pAlgorithm the address of this master algorithm
     */
    virtual void FindAmbiguousPfos(
        const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TAGGING_BASE_TOOL_H
