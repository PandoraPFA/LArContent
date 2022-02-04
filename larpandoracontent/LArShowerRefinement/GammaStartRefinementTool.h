/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
#define LAR_GAMMA_START_REFINEMENT_TOOL_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

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

    void BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType tpcView, ProtoShowerVector &protoShowerVector);


    void FillTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

    void AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const AngularDecompositionMap &angularDecompositionMapU, 
                             const AngularDecompositionMap &angularDecompositionMapV, const AngularDecompositionMap &angularDecompositionMapW);

    void ObtainAngularPeakVector(ShowerStartRefinementAlgorithm *const pAlgorithm, AngularDecompositionMap &angularDecompositionMapU, 
        AngularDecompositionMap &angularDecompositionMapV, AngularDecompositionMap &angularDecompositionMapW, AngularPeakVector &angularPeakVector);

    void FillMCParticleToHitMap(const pandora::CaloHitList *const pCaloHitList, const pandora::HitType tpcView, LArMCParticleHelper::MCContributionMap &mcParticleToHitMap);

    void FillOutPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, pandora::CaloHitList &unavailableHits, ProtoShowerVector &protoShowerVector);

    void RemoveTrackPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVector);

    bool IsTrack(const ProtoShower &protoShower);

    void RemoveConnectionPathway(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveConnectionPathway(const ProtoShower &protoShower);

    int m_counter;
    int m_showerCounter;
    float m_trackSearchWindow;
    float m_minTrackBlipMean;
    float m_electronFraction;
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
