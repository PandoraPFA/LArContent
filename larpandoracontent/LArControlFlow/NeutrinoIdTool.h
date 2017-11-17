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
class NeutrinoIdTool : public NeutrinoIdBaseTool
{
public:
    void SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_ID_TOOL_H
