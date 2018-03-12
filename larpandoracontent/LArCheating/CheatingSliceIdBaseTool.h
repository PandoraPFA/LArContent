/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceIdBaseTool.h
 *
 *  @brief  Header file for the cheating slice id base tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SLICE_ID_BASE_TOOL_H
#define LAR_CHEATING_SLICE_ID_BASE_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include <functional>

namespace lar_content
{

/**
 *  @brief  CheatingSliceIdBaseTool class
 */
class CheatingSliceIdBaseTool : public SliceIdBaseTool
{
public:
    virtual void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos) = 0;

    /**
     *  @brief  Get the target particle weight in a list of pfos
     *
     *  @param  pPfoList address of the pfo list
     *  @param  objectOwnedByMaster whether the pfos (and hits) are owned by the master Pandora instance, or a worker instance
     *  @param  targetParticleWeight the target particle weight
     *  @param  totalWeight the total weight
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     */
    static void GetTargetParticleWeight(const pandora::PfoList *const pPfoList, const bool objectOwnedByMaster, float &targetParticleWeight, float &totalWeight, std::function<bool(const pandora::MCParticle *const)> fCriteria);

    /**
     *  @brief  Get the target particle weight for a calo hit
     *
     *  @param  pCaloHit address of the calo hit
     *  @param  objectOwnedByMaster whether the hit is owned by the master Pandora instance, or a worker instance
     *  @param  targetParticleWeight the target particle weight
     *  @param  totalWeight the total weight
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     */
    static void GetTargetParticleWeight(const pandora::CaloHit *const pCaloHit, const bool objectOwnedByMaster, float &targetParticleWeight, float &totalWeight, std::function<bool(const pandora::MCParticle *const)> fCriteria);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SLICE_ID_BASE_TOOL_H
