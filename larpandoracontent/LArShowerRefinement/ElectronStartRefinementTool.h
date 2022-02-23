/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_ELECTRON_START_REFINEMENT_TOOL_H
#define LAR_ELECTRON_START_REFINEMENT_TOOL_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

namespace lar_content
{

/**
 *  @brief  ElectronStartRefinementTool class
 */
class ElectronStartRefinementTool : public ShowerStartRefinementBaseTool
{
public:
    ElectronStartRefinementTool();
    ~ElectronStartRefinementTool();

    bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition);

    void BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType hitType, ElectronProtoShowerVector &protoShowerVector,
        pandora::CaloHitList &usedCaloHitList);

    void BuildHelperProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, 
        const pandora::HitType tpcView, pandora::CartesianPointVector &significantPeakDirections, pandora::CaloHitList &usedHitList);

    void CollectHitsWithinROI(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerHitList, 
        const pandora::CaloHitList &allHitList, const pandora::CartesianVector &projectedNuVertexPosition, pandora::CaloHitList &collectedHits);

    void GetAngularExtrema(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &viewHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, float &lowestTheta, float &highestTheta);

    void CollectHitsWithinExtrema(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &viewHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, const float lowestTheta, const float highestTheta, pandora::CaloHitList &collectedHits);

    pandora::CartesianVector GetShowerVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const pandora::CartesianVector &projectedNuVertexPosition);

    TwoDSlidingShowerFitResult FitShower(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType);

    bool IsSpineCoincident(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &projectedNuVertexPosition, 
        const pandora::CartesianVector &peakDirection, const pandora::CaloHitList &viewShowerHitList, const pandora::CartesianVector &showerVertexPosition, 
        const pandora::CaloHitList &showerSpineHitList);

    bool IsShowerExtendable(ShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
        const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW);

    bool IsElectronPathway(ShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
        const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW);

    bool ArePathwaysConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW);

   bool AreShowerStartsConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float allowance);

    bool AreDirectionsConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float allowance);

    void ExtendShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const ElectronProtoShowerVector &protoShowerVectorU, const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW);

    void ExtendShowerOneView(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShowerVector &showersToExtend);

    void ExtendShowerTwoView(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShowerVector &showersToExtend);

    void ExtendShowerThreeView(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShowerVector &showersToExtend);

    void RefineHitsToAdd(ShowerStartRefinementAlgorithm *const pAlgorithm, ElectronProtoShower &protoShower,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType hitType, ElectronProtoShowerVector &protoShowerVector, 
        const pandora::CartesianPointVector &significantPeakDirections);

    pandora::CaloHitList FindContinuousPath(const pandora::CaloHitList &refinedHitList, const pandora::CartesianVector &nuVertexPosition);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_showerSlidingFitWindow;
    float m_maxCoincideneTransverseSeparation;
    float m_minSpinePurity;
    float m_maxAngularDeviation;
    float m_maxXSeparation;
    bool m_extendMode;
    bool m_moveVertexMode;
};

} // namespace lar_content

#endif // #ifndef LAR_ELECTRON_START_REFINEMENT_TOOL_H
