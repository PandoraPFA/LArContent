/**
 *  @file   larpandoracontent/LArCheating/CheatingBeamParticleIdTool.h
 *
 *  @brief  Header file for the cheating beam particle id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_BEAM_PARTICLE_ID_TOOL_H
#define LAR_CHEATING_BEAM_PARTICLE_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingBeamParticleIdTool class
 */
class CheatingBeamParticleIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Constructor
     */
    CheatingBeamParticleIdTool();

    /**
     *  @brief  Get the beam particle weight for a list of pfos
     *
     *  @param  pPfoList address of the pfo list
     *  @param  objectOwnedByMaster whether the pfos (and hits) are owned by the master Pandora instance, or a worker instance
     *  @param  beamParticleWeight the beam particle weight
     *  @param  totalWeight the total weight
     */
    static void GetBeamParticleWeight(const pandora::PfoList *const pPfoList, const bool objectOwnedByMaster, float &beamParticleWeight, float &totalWeight);

    /**
     *  @brief  Get the beam particle weight for a calo hit
     *
     *  @param  pCaloHit address of the calo hit
     *  @param  objectOwnedByMaster whether the hit is owned by the master Pandora instance, or a worker instance
     *  @param  beamParticleWeight the beam particle weight
     *  @param  totalWeight the total weight
     */
    static void GetBeamParticleWeight(const pandora::CaloHit *const pCaloHit, const bool objectOwnedByMaster, float &beamParticleWeight, float &totalWeight);

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    float     m_minWeightFraction;     ///< The minimum weight fraction for identifying a slice as a beam particle 
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_BEAM_PARTICLE_ID_TOOL_H
