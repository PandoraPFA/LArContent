/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoIdTool.h
 *
 *  @brief  Header file for the neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_ID_TOOL_H
#define LAR_NEUTRINO_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoIdTool class
 */
class NeutrinoIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoIdTool();

    void SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool        m_selectAllNeutrinos;               ///< First approach: select all neutrinos, as opposed to selecting all cosmics
    bool        m_selectOnlyFirstSliceNeutrinos;    ///< First approach: select first slice neutrinos, cosmics for all subsequent slices
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_ID_TOOL_H
