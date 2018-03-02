/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoIdTool.h
 *
 *  @brief  Header file for the cheating neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_ID_TOOL_H
#define LAR_CHEATING_NEUTRINO_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoIdTool class
 */
class CheatingNeutrinoIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Get the neutrino weight for a list of pfos
     *
     *  @param  pPfoList address of the pfo list
     *  @param  objectOwnedByMaster whether the pfos (and hits) are owned by the master Pandora instance, or a worker instance
     *  @param  neutrinoWeight the neutrino weight
     *  @param  totalWeight the total weight
     */
    static void GetNeutrinoWeight(const pandora::PfoList *const pPfoList, const bool objectOwnedByMaster, float &neutrinoWeight, float &totalWeight);

    /**
     *  @brief  Get the neutrino weight for a calo hit
     *
     *  @param  pCaloHit address of the calo hit
     *  @param  objectOwnedByMaster whether the hit is owned by the master Pandora instance, or a worker instance
     *  @param  neutrinoWeight the neutrino weight
     *  @param  totalWeight the total weight
     */
    static void GetNeutrinoWeight(const pandora::CaloHit *const pCaloHit, const bool objectOwnedByMaster, float &neutrinoWeight, float &totalWeight);

    void SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_ID_TOOL_H
