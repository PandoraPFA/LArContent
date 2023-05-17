/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h
 *
 *  @brief  Header file for the peak direction finder tool class.
  *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_FINDER_TOOL_H
#define LAR_SHOWER_START_FINDER_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

namespace lar_content
{

class ShowerStartFinderTool : public pandora::AlgorithmTool
{
public:
    ShowerStartFinderTool();

    pandora::StatusCode Run(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &peakDirection, const pandora::HitType hitType, 
        const pandora::CaloHitList &showerSpineHitList, pandora::CartesianVector &showerStartPosition, pandora::CartesianVector &showerStartDirection);

private:
    typedef std::map<const pandora::CaloHit*, float> LongitudinalPositionMap;
    typedef std::map<int, float> EnergySpectrumMap;
    typedef std::map<int, pandora::CaloHitList> LayerToHitMap;

    void ObtainLongitudinalDecomposition(const TwoDSlidingFitResult &spineTwoDSlidingFit, const pandora::CaloHitList &showerSpineHitList, 
        LongitudinalPositionMap &longitudinalPositionMap) const;

    void GetEnergyDistribution(const pandora::CaloHitList &showerSpineHitList, const LongitudinalPositionMap &longitudinalPositionMap, 
        EnergySpectrumMap &energySpectrumMap) const;

    pandora::StatusCode FindShowerStartAndDirection(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const TwoDSlidingFitResult &spineTwoDSlidingFit, const EnergySpectrumMap &energySpectrumMap, const pandora::CaloHitList &showerSpineHitList, 
        const bool isEndDownstream, pandora::CartesianVector &showerStartPosition, pandora::CartesianVector &showerStartDirection) const;

    void CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, 
        float &energySigma) const;

    bool IsShowerTopology(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, const TwoDSlidingFitResult &spineTwoDSlidingFit, 
        const float longitudinalDistance, const pandora::CaloHitList &showerSpineHitList, const bool isEndDownstream) const;

    void ConvertLongitudinalProjectionToGlobal(const TwoDSlidingFitResult &spineTwoDSlidingFit, const float longitudinalDistance, 
        pandora::CartesianVector &globalPosition, pandora::CartesianVector &globalDirection) const;

    pandora::StatusCode BuildShowerRegion(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, const pandora::CaloHitList &showerSpineHitList, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianPointVector &showerRegionPositionVector) const;

    pandora::StatusCode CharacteriseShowerTopology(const pandora::CartesianPointVector &showerRegionPositionVector, const pandora::CartesianVector &showerStartPosition,
        const bool isEndDownstream, const pandora::CartesianVector &showerStartDirection, pandora::CartesianVector &positiveEdgeStart, pandora::CartesianVector &positiveEdgeEnd, 
        pandora::CartesianVector &negativeEdgeStart, pandora::CartesianVector &negativeEdgeEnd, bool &isBetween, bool &doesStraddle) const;

    bool IsClockwiseRotation(const pandora::CartesianVector &showerStartDirection, const pandora::CartesianVector &displacementVector) const;

    pandora::StatusCode GetBoundaryExtremalPoints(const TwoDSlidingShowerFitResult &showerTwoDSlidingFit, const LayerFitResultMap &layerFitResultMap, 
        const pandora::CartesianVector &showerStartPosition, const int showerStartLayer, const int showerEndLayer, pandora::CartesianVector &boundaryEdgeStart, 
        pandora::CartesianVector &boundaryEdgeEnd) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_spineSlidingFitWindow;
    float m_longitudinalCoordinateBinSize;
    unsigned int m_nInitialEnergyBins;
    float m_minSigmaDeviation;
    float m_minShowerOpeningAngle;
    float m_molliereRadius;
    unsigned int m_showerSlidingFitWindow;
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_FINDER_TOOL_H
