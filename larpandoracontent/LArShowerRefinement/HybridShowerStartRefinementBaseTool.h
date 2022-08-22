/**
 *  @file   larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementBaseTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_HYBRID_SHOWER_START_REFINEMENT_BASE_TOOL_H
#define LAR_HYBRID_SHOWER_START_REFINEMENT_BASE_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementAlgorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

class HybridShowerStartRefinementBaseTool : public pandora::AlgorithmTool
{
public:
    HybridShowerStartRefinementBaseTool();

    virtual bool Run(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition, 
        const pandora::CaloHitList *const pCaloHitListU, const pandora::CaloHitList *const pCaloHitListV, const pandora::CaloHitList *const pCaloHitListW) = 0;

    typedef std::map<int, float> AngularDecompositionMap;
    typedef std::map<const pandora::CaloHit*, float> LongitudinalPositionMap;
    typedef std::map<int, float> EnergySpectrumMap;
    typedef std::map<int, pandora::CaloHitList> LayerToHitMap;

protected:
    bool HasPathToNuVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &neutrinoVertex) const;

    bool HasPathToNuVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &neutrinoVertex, const pandora::HitType &hitType) const;

    void FillAngularDecompositionMap(const pandora::CaloHitList &viewShowerHitList, const pandora::CartesianVector &projectedNuVertexPosition, 
        AngularDecompositionMap &angularDecompositionMap);

    void SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap);

    void ObtainPeakVector(AngularDecompositionMap &angularDecompositionMap, pandora::IntVector &viewPeakVector);

    bool FindBestAngularPeak(AngularDecompositionMap &angularDecompositionMap, pandora::IntVector &viewPeakVector, pandora::IntVector &investigatedPeaks, int &bestAngularPeak);

    void FindShowerSpine(const HybridShowerStartRefinementAlgorithm *pAlgorithm, const pandora::CaloHitList &viewShowerHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, const pandora::CartesianVector &initialDirection, pandora::CaloHitList &unavailableHitList, 
        pandora::CaloHitList &showerSpineHitList);

    bool CollectSubsectionHits(const HybridShowerStartRefinementAlgorithm *pAlgorithm, const TwoDSlidingFitResult &extrapolatedFit, const pandora::CartesianVector &extrapolatedStartPosition, 
       const pandora::CartesianVector &extrapolatedEndPosition, const pandora::CartesianVector &extrapolatedDirection, const bool isEndDownstream, const pandora::CaloHitList &viewShowerHitList, 
       pandora::CartesianPointVector &runningFitPositionVector, pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList);

    void CollectConnectedHits(const HybridShowerStartRefinementAlgorithm *pAlgorithm, const pandora::CaloHitList &collectedHits, const pandora::CartesianVector &extrapolatedStartPosition, 
        const pandora::CartesianVector &extrapolatedDirection, pandora::CartesianPointVector &runningFitPositionVector, pandora::CaloHitList &unavailableHitList, 
        pandora::CaloHitList &showerSpineHitList);

    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineDirection, 
        const float distanceToLine) const;

    float GetClosestDistance(const pandora::CartesianVector &position, const pandora::CartesianPointVector &testPositions) const;

    void ObtainLongitudinalDecomposition(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        LongitudinalPositionMap &longitudinalPositionMap);

    void GetEnergyDistribution(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        const LongitudinalPositionMap &longitudinalPositionMap, EnergySpectrumMap &energySpectrumMap);

    bool IsShowerTopology(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, const pandora::CaloHitList &showerPfoHits, 
        const pandora::CaloHitList &showerSpineHits, const bool isEndDownstream);

    void CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, float &energySigma);

    bool FindShowerStart(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &projectedNuVertexPosition, const pandora::CartesianVector &peakDirection, 
        const LongitudinalPositionMap &longitudinalPositionMap, const EnergySpectrumMap &energySpectrumMap, const pandora::CaloHitList &showerSpineHitList, 
        pandora::CartesianVector &showerStartPosition, const pandora::CaloHitList &showerPfoHitList, const bool isEndDownstream, ProtoShowerVector &protoShowerVector,
        const bool isHelper);

    void ConvertLongitudinalProjectionToGlobalPosition(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, const float longitudinalDistance, 
        pandora::CartesianVector &globalPosition, pandora::CartesianVector &globalDirection);

    pandora::StatusCode FillHaloHitPositionVector(const pandora::CaloHitList &viewShowerHitList, const pandora::CaloHitList &showerSpineHitList, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianPointVector &haloHitPositionVector);

    pandora::StatusCode CharacteriseShower(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, const pandora::CaloHitList &showerSpineHits, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianVector &positiveEdgeStart, pandora::CartesianVector &positiveEdgeEnd, pandora::CartesianVector &negativeEdgeStart, pandora::CartesianVector &negativeEdgeEnd, 
        bool &isBetween, bool &doesStraddle);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxDistanceForConnection;
    float m_pathwaySearchRegion;
    float m_theta0XZBinSize;
    int m_smoothingWindow;
    float m_growingFitInitialLength;
    float m_macroSlidingFitWindow;
    float m_growingFitSegmentLength;
    float m_distanceToLine;
    float m_initialFitDistanceToLine;
    int m_maxFittingHits;
    float m_longitudinalCoordinateBinSize;
    float m_hitConnectionDistance;
    unsigned int m_minInitialHitsFound;
    int m_microSlidingFitWindow;
    unsigned int m_nInitialEnergyBins;
    float m_minSigmaDeviation;
    float m_molliereRadius;
    int m_showerSlidingFitWindow;
    float m_minShowerOpeningAngle;

    //private:
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
