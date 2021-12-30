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

    void FillAngularDecompositionMap(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &viewShowerHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, AngularDecompositionMap &angularDecompositionMap);

    void SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap);

    void ObtainPeakVector(AngularDecompositionMap &angularDecompositionMap, pandora::IntVector &viewPeakVector);

    void GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        const LongitudinalPositionMap &longitudinalPositionMap, EnergySpectrumMap &energySpectrumMap);

    void ObtainLongitudinalDecomposition(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        LongitudinalPositionMap &longitudinalPositionMap);

    void CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, float &energySigma);

    bool FindShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &projectedNuVertexPosition, const pandora::CartesianVector &peakDirection, 
        const LongitudinalPositionMap &longitudinalPositionMap, const EnergySpectrumMap &energySpectrumMap, const pandora::CaloHitList &showerSpineHitList, 
        pandora::CartesianVector &showerStartPosition, const pandora::CaloHitList &showerPfoHitList, const bool isEndDownstream, ProtoShowerVector &protoShowerVector);

    void ConvertLongitudinalProjectionToGlobalPosition(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, const float longitudinalDistance, 
        pandora::CartesianVector &globalPosition, pandora::CartesianVector &globalDirection);

    pandora::StatusCode FillHaloHitPositionVector(const pandora::CaloHitList &viewShowerHitList, const pandora::CaloHitList &showerSpineHitList, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianPointVector &haloHitPositionVector);


    void FillTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

    void AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const AngularDecompositionMap &angularDecompositionMapU, 
                             const AngularDecompositionMap &angularDecompositionMapV, const AngularDecompositionMap &angularDecompositionMapW);

    void ObtainAngularPeakVector(ShowerStartRefinementAlgorithm *const pAlgorithm, AngularDecompositionMap &angularDecompositionMapU, 
        AngularDecompositionMap &angularDecompositionMapV, AngularDecompositionMap &angularDecompositionMapW, AngularPeakVector &angularPeakVector);

    bool FindBestAngularPeak(AngularDecompositionMap &angularDecompositionMap, pandora::IntVector &viewPeakVector, pandora::IntVector &investigatedPeaks, int &bestAngularPeak);


    bool IsShowerTopology(ShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, const pandora::CaloHitList &showerPfoHits, 
        const pandora::CaloHitList &showerSpineHits, const bool isEndDownstream);

    pandora::StatusCode CharacteriseShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, const pandora::CaloHitList &showerSpineHits, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianVector &positiveEdgeStart, pandora::CartesianVector &positiveEdgeEnd, pandora::CartesianVector &negativeEdgeStart, pandora::CartesianVector &negativeEdgeEnd, 
        bool &isBetween, bool &doesStraddle);

    void FillMCParticleToHitMap(const pandora::CaloHitList *const pCaloHitList, const pandora::HitType tpcView, LArMCParticleHelper::MCContributionMap &mcParticleToHitMap);

    void FillOutPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, pandora::CaloHitList &unavailableHits, ProtoShowerVector &protoShowerVector);

    void RemoveTrackPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVector);

    bool IsTrack(const ProtoShower &protoShower);

    void RemoveConnectionPathway(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveConnectionPathway(const ProtoShower &protoShower);

    int m_counter;
    float m_theta0XZBinSize;
    float m_pathwaySearchRegion;
    int m_smoothingWindow;
    int m_showerCounter;
    int m_microSlidingFitWindow;
    float m_minSigmaDeviation;
    float m_trackSearchWindow;
    unsigned int m_nInitialEnergyBins;
    float m_minTrackBlipMean;
    int m_showerSlidingFitWindow;
    float m_molliereRadius;
    float m_minShowerOpeningAngle;
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
