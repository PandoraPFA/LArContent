/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoSliceSelectionTool.h
 *
 *  @brief  Header file for the neutrino slice selection tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_SLICE_SELECTION_TOOL_H
#define LAR_NEUTRINO_SLICE_SELECTION_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoSliceSelectionTool class
 */
class NeutrinoSliceSelectionTool : public SliceSelectionBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoSliceSelectionTool();

    void SelectSlices(const pandora::Algorithm *const pAlgorithm, const SliceVector &inputSliceVector,
        SliceVector &outputSliceVector);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    int             m_maxSlices;                ///< The maximum number of slices to retain (0 to retain all) - default 0
    float           m_threshold;                ///< The minimum cut threshold to retain a slice (< 0 for no threshold) - default -1
    std::string     m_cutVariable;              ///< The variable to cut on ("purity" or "completeness") - default "completeness"
    std::string     m_mcParticleListName;       ///< The name for the MCParticleList associated with the slice selection
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_SLICE_SELECTION_TOOL_H

