/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
#define LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

class ShowerCore
{
public:
    ShowerCore();

    ShowerCore(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection, const pandora::CaloHitList &coreHitList);

    void SetStartPosition(const pandora::CartesianVector &startPosition);

    void SetStartDirection(const pandora::CartesianVector &startDirection);

    // when you're not feeling lazy, change this to private 
    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
    pandora::CaloHitList m_coreHitList;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection, const pandora::CaloHitList &coreHitList) : 
    m_startPosition(startPosition), m_startDirection(startDirection), m_coreHitList(coreHitList)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f), m_coreHitList()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerCore::SetStartPosition(const pandora::CartesianVector &startPosition)
{
    m_startPosition = startPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerCore::SetStartDirection(const pandora::CartesianVector &startDirection)
{
    m_startDirection = startDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ConnectionPathway
{
public:
    ConnectionPathway();

    ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection, const pandora::CaloHitList &pathwayHitList);

    // when you're not feeling lazy, change this to private
    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
    pandora::CaloHitList m_pathwayHitList;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection, const pandora::CaloHitList &pathwayHitList) : 
    m_startPosition(startPosition), m_startDirection(startDirection), m_pathwayHitList(pathwayHitList)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f), m_pathwayHitList()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ProtoShower
{
public: 
    ProtoShower();

    ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const bool isHelper);

    // when you're not feeling lazy, change this to private  
    ShowerCore m_showerCore;
    ConnectionPathway m_connectionPathway;
    bool m_isHelper;
};

typedef std::vector<ProtoShower> ProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const bool isHelper) : 
    m_showerCore(showerCore), m_connectionPathway(connectionPathway), m_isHelper(isHelper)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower() : m_showerCore(), m_connectionPathway(), m_isHelper(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ElectronProtoShower : public ProtoShower
{
public: 
    ElectronProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const bool isHelper, const pandora::CaloHitList &hitsToAdd);

    // when you're not feeling lazy, change this to private  
    pandora::CaloHitList m_hitsToAdd;
};

typedef std::vector<ElectronProtoShower> ElectronProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ElectronProtoShower::ElectronProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const bool isHelper, 
    const pandora::CaloHitList &hitsToAdd) : ProtoShower(showerCore, connectionPathway, isHelper), m_hitsToAdd(hitsToAdd)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class AngularPeak
{
public:
    AngularPeak(int uPeakBin, int vPeakBin, int wPeakBin, float chiSquared, float binWeightSum);

    int m_uPeakBin;
    int m_vPeakBin;
    int m_wPeakBin;
    float m_chiSquared;
    float m_binWeightSum;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline AngularPeak::AngularPeak(int uPeakBin, int vPeakBin, int wPeakBin, float chiSquared, float binWeightSum)
{
    m_uPeakBin = uPeakBin;
    m_vPeakBin = vPeakBin;
    m_wPeakBin = wPeakBin;
    m_chiSquared = chiSquared;
    m_binWeightSum = binWeightSum;
}

typedef std::vector<AngularPeak> AngularPeakVector;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ShowerStartRefinementBaseTool : public pandora::AlgorithmTool
{
public:
    ShowerStartRefinementBaseTool();

    virtual bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition) = 0;

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

    void FindShowerSpine(const ShowerStartRefinementAlgorithm *pAlgorithm, const pandora::CaloHitList &viewShowerHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition, const pandora::CartesianVector &initialDirection, pandora::CaloHitList &unavailableHitList, 
        pandora::CaloHitList &showerSpineHitList);

    bool CollectSubsectionHits(const ShowerStartRefinementAlgorithm *pAlgorithm, const TwoDSlidingFitResult &extrapolatedFit, const pandora::CartesianVector &extrapolatedStartPosition, 
       const pandora::CartesianVector &extrapolatedEndPosition, const pandora::CartesianVector &extrapolatedDirection, const bool isEndDownstream, const pandora::CaloHitList &viewShowerHitList, 
       pandora::CartesianPointVector &runningFitPositionVector, pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList);

    void CollectConnectedHits(const ShowerStartRefinementAlgorithm *pAlgorithm, const pandora::CaloHitList &collectedHits, const pandora::CartesianVector &extrapolatedStartPosition, 
        const pandora::CartesianVector &extrapolatedDirection, pandora::CartesianPointVector &runningFitPositionVector, pandora::CaloHitList &unavailableHitList, 
        pandora::CaloHitList &showerSpineHitList);

    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineDirection, 
        const float distanceToLine) const;

    float GetClosestDistance(const pandora::CartesianVector &position, const pandora::CartesianPointVector &testPositions) const;

    void ObtainLongitudinalDecomposition(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        LongitudinalPositionMap &longitudinalPositionMap);

    void GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, 
        const LongitudinalPositionMap &longitudinalPositionMap, EnergySpectrumMap &energySpectrumMap);

    bool IsShowerTopology(ShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, const pandora::CaloHitList &showerPfoHits, 
        const pandora::CaloHitList &showerSpineHits, const bool isEndDownstream);

    void CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, float &energySigma);

    bool FindShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &projectedNuVertexPosition, const pandora::CartesianVector &peakDirection, 
        const LongitudinalPositionMap &longitudinalPositionMap, const EnergySpectrumMap &energySpectrumMap, const pandora::CaloHitList &showerSpineHitList, 
        pandora::CartesianVector &showerStartPosition, const pandora::CaloHitList &showerPfoHitList, const bool isEndDownstream, ProtoShowerVector &protoShowerVector,
        const bool isHelper);

    void ConvertLongitudinalProjectionToGlobalPosition(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerSpineHitList, const float longitudinalDistance, 
        pandora::CartesianVector &globalPosition, pandora::CartesianVector &globalDirection);

    pandora::StatusCode FillHaloHitPositionVector(const pandora::CaloHitList &viewShowerHitList, const pandora::CaloHitList &showerSpineHitList, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianPointVector &haloHitPositionVector);

    pandora::StatusCode CharacteriseShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, const pandora::CaloHitList &showerSpineHits, 
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

        class SortByDistanceToPoint
        {
        public:
            /**
             *  @brief  Constructor
             *
             *  @param  referencePoint the point relative to which constituent hits are ordered
             */
        SortByDistanceToPoint(const pandora::CartesianVector referencePoint) : m_referencePoint(referencePoint)
            {
            }

            /**
             *  @brief  Sort constituent hits by their position relative to a referencePoint
             *
             *  @param  lhs first constituent hit
             *  @param  rhs second constituent hit
             *
             *  @return  whether lhs hit is closer to the referencePoint than the rhs hit
             */
            bool operator()(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);
            bool operator()(const pandora::CaloHit *const lhs, const pandora::CaloHit *const rhs);

        private:
            const pandora::CartesianVector m_referencePoint; ///< The point relative to which constituent hits are ordered
        };

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
