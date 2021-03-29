/**
 *  @file   larpandoracontent/LArCheating/CheatingBeamParticleIdTool.h
 *
 *  @brief  Header file for the cheating beam particle id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_BEAM_PARTICLE_ID_TOOL_H
#define LAR_CHEATING_BEAM_PARTICLE_ID_TOOL_H 1

#include "larpandoracontent/LArCheating/CheatingSliceIdBaseTool.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingBeamParticleIdTool class
 */
class CheatingBeamParticleIdTool : public CheatingSliceIdBaseTool
{
public:
    /**
     *  @brief  Constructor
     */
    CheatingBeamParticleIdTool();

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minWeightFraction; ///< The minimum weight fraction for identifying a slice as a beam particle
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_BEAM_PARTICLE_ID_TOOL_H
