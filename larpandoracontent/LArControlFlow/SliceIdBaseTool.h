/**
 *  @file   larpandoracontent/LArControlFlow/SliceIdBaseTool.h
 *
 *  @brief  Header file for the stitching tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_ID_BASE_TOOL_H
#define LAR_SLICE_ID_BASE_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

typedef std::vector<pandora::CaloHitList> SliceVector;
typedef std::vector<pandora::PfoList> SliceHypotheses;
/**
 *  @brief  SliceIdBaseTool class
 */
class SliceIdBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Select which reconstruction hypotheses to use; neutrino outcomes or cosmic-ray muon outcomes for each slice
     *
     *  @param  pAlgorithm the address of the master instance, used to access MCParticles when in training mode
     *  @param  nuSliceHypotheses the parent pfos representing the neutrino outcome for each slice
     *  @param  crSliceHypotheses the parent pfos representing the cosmic-ray muon outcome for each slice
     *  @param  sliceNuPfos to receive the list of selected pfos
     */
    virtual void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos) = 0;
};

}

#endif // #ifndef LAR_SLICE_ID_BASE_TOOL_H
