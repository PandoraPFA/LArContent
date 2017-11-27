/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleIdTool.h
 *
 *  @brief  Header file for the beam particle id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_BEAM_PARTICLE_ID_TOOL_H
#define LAR_BEAM_PARTICLE_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  BeamParticleIdTool class
 */
class BeamParticleIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    BeamParticleIdTool();

    void SelectOutputPfos(const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool        m_selectAllBeamParticles;               ///< First approach: select all beam particles, as opposed to selecting all cosmics
    bool        m_selectOnlyFirstSliceBeamParticles;    ///< First approach: select first slice beam particles, cosmics for all subsequent slices
};

} // namespace lar_content

#endif // #ifndef LAR_BEAM_PARTICLE_ID_TOOL_H
