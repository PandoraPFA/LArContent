/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
#define LAR_GAMMA_START_REFINEMENT_TOOL_H 1

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

    void FillTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

    void FillAngularDecompositionMap(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType hitType, DeviationAngleMap &deviationAngleMap);

    void SmoothAngularDecompositionMap(DeviationAngleMap &deviationAngleMap);

    void ObtainAngularPeakVector(ShowerStartRefinementAlgorithm *const pAlgorithm, DeviationAngleMap &deviationAngleMapU, 
        DeviationAngleMap &deviationAngleMapV, DeviationAngleMap &deviationAngleMapW, AngularPeakVector &angularPeakVector);

    void ObtainViewPeakVector(DeviationAngleMap &deviationAngleMap, pandora::IntVector &viewPeakVector);

    int FindBestAngularPeak(DeviationAngleMap &deviationAngleMap, pandora::IntVector &viewPeakVector);

    void GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, const LongitudinalPositionMap &longitudinalPositionMap);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveConnectionPathway(const ProtoShower &protoShower);

    int m_counter;
    int m_showerCounter;
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
