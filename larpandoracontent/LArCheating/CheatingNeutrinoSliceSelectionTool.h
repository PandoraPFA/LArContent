/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoSliceSelectionTool.h
 *
 *  @brief  Header file for the neutrino slice selection tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_SLICE_SELECTION_TOOL_H
#define LAR_CHEATING_NEUTRINO_SLICE_SELECTION_TOOL_H 1

#include "larpandoracontent/LArCheating/CheatingSliceSelectionTool.h"

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoSliceSelectionTool class
 */
class CheatingNeutrinoSliceSelectionTool : public CheatingSliceSelectionTool
{
protected:
    /**
     *  @brief  Template method to determine if an MC particle matches the target criteria for slice selection. Return true if match.
     *
     *  @param  mcParticle the MC particle to check
     */
    bool IsTarget(const pandora::MCParticle *const mcParticle) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_SLICE_SELECTION_TOOL_H
