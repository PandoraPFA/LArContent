/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
#define LAR_GAMMA_START_REFINEMENT_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

namespace lar_content
{

/**
 *  @brief  GammaStartRefinementTool class
 */
class GammaStartRefinementTool : public ShowerStartRefinementBaseTool
{
public:
    GammaStartRefinementTool();
    ~GammaStartRefinementTool();

    bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveConnectionPathway(const ProtoShower &protoShower);

    int m_counter;
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
