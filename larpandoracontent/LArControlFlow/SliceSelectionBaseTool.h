/**
 *  @file   larpandoracontent/LArControlFlow/SliceSelectionBaseTool.h
 *
 *  @brief  Header file for the stitching tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_SELECTION_BASE_TOOL_H
#define LAR_SLICE_SELECTION_BASE_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

typedef std::vector<pandora::CaloHitList> SliceVector;
typedef std::vector<pandora::PfoList> SliceHypotheses;
/**
 *  @brief  SliceSelectionBaseTool class
 */
class SliceSelectionBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Select which slice(s) to use; neutrino or beam slices
     *
     *  @param  pAlgorithm the address of the master instance, used to access MCParticles when in training mode
     *  @param  inputSliceVector the initial slice vector
     *  @param  outputSliceVector the output slice vector
     */
    virtual void SelectSlices(const pandora::Algorithm *const pAlgorithm, const SliceVector &inputSliceVector, SliceVector &outputSliceVector) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_SLICE_SELECTION_BASE_TOOL_H
