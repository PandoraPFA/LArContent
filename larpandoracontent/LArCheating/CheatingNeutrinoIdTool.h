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
    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    /**
     *  @brief  Get the neutrino weight for a CaloHit
     *
     *  @param  pCaloHit address of the CaloHit
     *  @param  neutrinoWeight the neutrino weight
     *  @param  totalWeight the total weight
     */
    void GetNeutrinoWeight(const pandora::CaloHit *const pCaloHit, float &neutrinoWeight, float &totalWeight) const;

    /**
     *  @brief  Get the neutrino weight for a list of PFOs
     *
     *  @param  pPfoList address of the PFO list
     *  @param  neutrinoWeight the neutrino weight
     *  @param  totalWeight the total weight
     */
    void GetNeutrinoWeight(const pandora::PfoList *const pPfoList, float &neutrinoWeight, float &totalWeight) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_ID_TOOL_H
