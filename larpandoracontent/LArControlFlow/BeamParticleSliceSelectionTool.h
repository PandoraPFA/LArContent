/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleSliceSelectionTool.h
 *
 *  @brief  Header file for the test beam slice selection tool class.
 *
 *  $Log: $
 */
#ifndef LAR_BEAM_PARTICLE_SLICE_SELECTION_TOOL_H
#define LAR_BEAM_PARTICLE_SLICE_SELECTION_TOOL_H 1

#include "larpandoracontent/LArControlFlow/GenericSliceSelectionTool.h"

namespace lar_content
{

/**
 *  @brief  BeamParticleSliceSelectionTool class
 */
class BeamParticleSliceSelectionTool : public GenericSliceSelectionTool
{
public:
    /**
     *  @brief  Default constructor
     */
    BeamParticleSliceSelectionTool();

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

#endif // #ifndef LAR_BEAM_PARTICLE_SLICE_SELECTION_TOOL_H

