/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool.
 *
 *  $Log: $
 */
#ifndef LAR_ELECTRON_START_REFINEMENT_TOOL_H
#define LAR_ELECTRON_START_REFINEMENT_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

namespace lar_content
{

/**
 *  @brief  ElectronStartRefinementTool class
 */
class ElectronStartRefinementTool : public ShowerStartRefinementBaseTool
{
public:
    ElectronStartRefinementTool();

    bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_ELECTRON_START_REFINEMENT_TOOL_H
