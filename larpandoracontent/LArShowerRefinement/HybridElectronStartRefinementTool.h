/**
 *  @file   larpandoracontent/LArShowerRefinement/HybridElectronStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_HYBRID_ELECTRON_START_REFINEMENT_TOOL_H
#define LAR_HYBRID_ELECTRON_START_REFINEMENT_TOOL_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"

#include "larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementBaseTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

namespace lar_content
{

/**
 *  @brief  HybridElectronStartRefinementTool class
 */
class HybridElectronStartRefinementTool : public HybridShowerStartRefinementBaseTool
{
public:
    HybridElectronStartRefinementTool();
    ~HybridElectronStartRefinementTool();

    bool Run(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition, const pandora::CaloHitList *const pCaloHitListU, const pandora::CaloHitList *const pCaloHitListV, 
        const pandora::CaloHitList *const pCaloHitListW);

    bool IsSensibleShower(const LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    void BuildProtoShowers(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType hitType, ElectronProtoShowerVector &protoShowerVector,
        pandora::CaloHitList &usedCaloHitList);

    void BuildHelperProtoShowers(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo ,const pandora::CartesianVector &nuVertexPosition, 
        const pandora::HitType tpcView, pandora::CartesianPointVector &significantPeakDirections, pandora::CaloHitList &usedHitList);

    void CollectHitsWithinROI(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerHitList, 
        const pandora::CaloHitList &allHitList, const pandora::CartesianVector &projectedNuVertexPosition, pandora::CaloHitList &collectedHits);

    void GetAngularExtrema(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &viewHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, float &lowestTheta, float &highestTheta);

    void CollectHitsWithinExtrema(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &viewHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, const float lowestTheta, const float highestTheta, pandora::CaloHitList &collectedHits);

    pandora::CartesianVector GetShowerVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const pandora::CartesianVector &projectedNuVertexPosition);

    TwoDSlidingShowerFitResult FitShower(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType);

    bool IsSpineCoincident(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &projectedNuVertexPosition, 
        const pandora::CartesianVector &peakDirection, const pandora::CaloHitList &viewShowerHitList, const pandora::CartesianVector &showerVertexPosition, 
        const pandora::CaloHitList &showerSpineHitList);

    bool IsShowerExtendable(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShower &protoShowerU, 
        const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW);

    bool IsElectronPathway(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
        const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW);

    bool ArePathwaysConsistent(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, 
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, Consistency &consistency);

    void ExtendShower(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW);

    void ExtendShowerOneView(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShowerVector &showersToExtend);

    void ExtendShowerTwoView(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShowerVector &showersToExtend);

    void ExtendShowerThreeView(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShowerVector &showersToExtend);

    void RefineHitsToAdd(HybridShowerStartRefinementAlgorithm *const pAlgorithm, ElectronProtoShower &protoShower,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType hitType, const pandora::CartesianPointVector &significantPeakDirections);

    bool IsShowerConnected(const pandora::CartesianVector &showerVertexPosition, const pandora::CartesianVector &projectedNuVertexPosition, 
        const pandora::CartesianVector &peakDirection);

    void AssignShowerHits(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::HitType &hitType, ElectronProtoShowerVector &protoShowerVector);

    pandora::CaloHitList FindContinuousPath(const pandora::CaloHitList &refinedHitList, const pandora::CartesianVector &nuVertexPosition);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_showerSlidingFitWindow;
    float m_maxCoincideneTransverseSeparation;
    float m_minSpinePurity;
    float m_maxAngularDeviation;
    float m_maxXSeparation;
    float m_maxSeparation;
    bool m_extendElectronMode;
    bool m_truncateGammaMode;
};

} // namespace lar_content

#endif // #ifndef LAR_ELECTRON_START_REFINEMENT_TOOL_H
