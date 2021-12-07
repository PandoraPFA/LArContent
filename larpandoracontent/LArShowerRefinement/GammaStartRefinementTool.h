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

    void AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const DeviationAngleMap &deviationAngleMapU, 
                             const DeviationAngleMap &deviationAngleMapV, const DeviationAngleMap &deviationAngleMapW);

    void FillAngularDecompositionMap(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType hitType, DeviationAngleMap &deviationAngleMap);

    void SmoothAngularDecompositionMap(DeviationAngleMap &deviationAngleMap);

    void ObtainAngularPeakVector(ShowerStartRefinementAlgorithm *const pAlgorithm, DeviationAngleMap &deviationAngleMapU, 
        DeviationAngleMap &deviationAngleMapV, DeviationAngleMap &deviationAngleMapW, AngularPeakVector &angularPeakVector);

    void ObtainViewPeakVector(DeviationAngleMap &deviationAngleMap, pandora::IntVector &viewPeakVector);

    bool FindBestAngularPeak(DeviationAngleMap &deviationAngleMap, pandora::IntVector &viewPeakVector, pandora::IntVector &investigatedPeaks, int &bestAngularPeak);

    void GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        pandora::IntVector &showerCounterVector, pandora::FloatVector &projectionVector, pandora::FloatVector &energies, EnergySpectrumMap &energySpectrumMap);

    void ObtainLongitudinalDecomposition(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        LongitudinalPositionMap &longitudinalPositionMap);

    bool FindShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const EnergySpectrumMap &energySpectrumMap, const pandora::CaloHitList &showerSpineHitList, pandora::CartesianVector &showerStartPosition, const pandora::CaloHitList &showerPfoHitList, const bool isEndDownstream);

    bool IsTrackBlip(const EnergySpectrumMap &energySpectrumMap, const EnergySpectrumMap::const_iterator &iter, const float initialMean, const float initialSigma);

    bool IsShowerTopology(ShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, const pandora::CaloHitList &showerPfoHits, const pandora::CaloHitList &showerSpineHits, const bool isEndDownstream);


    pandora::StatusCode CharacteriseShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, const pandora::CaloHitList &showerSpineHits, const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, float showerLengthFraction, pandora::CartesianVector &positiveEdgeStart, pandora::CartesianVector &positiveEdgeEnd, pandora::CartesianVector &positiveEdgeDirection, pandora::CartesianVector &negativeEdgeStart, pandora::CartesianVector &negativeEdgeEnd, pandora::CartesianVector &negativeEdgeDirection, bool &isBetween, bool &doesStraddle);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveConnectionPathway(const ProtoShower &protoShower);

    int m_counter;
    int m_showerCounter;
    int m_microSlidingFitWindow;
    float m_minSigmaDeviation;
    float m_trackSearchWindow;
    int m_nInitialEnergyBins;
    float m_minTrackBlipMean;
    int m_showerSlidingFitWindow;
    float m_molliereRadius;
    float m_minShowerOpeningAngle;
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
