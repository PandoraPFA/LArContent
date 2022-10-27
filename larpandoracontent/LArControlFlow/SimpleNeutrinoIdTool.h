/**
 *  @file   larpandoracontent/LArControlFlow/SimpleNeutrinoIdTool.h
 *
 *  @brief  Header file for the neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_NEUTRINO_ID_TOOL_H
#define LAR_SIMPLE_NEUTRINO_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/SliceIdBaseTool.h"

namespace lar_content
{

/**
 *  @brief  SimpleNeutrinoIdTool class
 */
class SimpleNeutrinoIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleNeutrinoIdTool();

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_selectAllNeutrinos;            ///< First approach: select all neutrinos, as opposed to selecting all cosmics
    bool m_selectOnlyFirstSliceNeutrinos; ///< First approach: select first slice neutrinos, cosmics for all subsequent slices
};

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_NEUTRINO_ID_TOOL_H
