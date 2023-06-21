/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h
 *
 *  @brief  Header file for the shower start finder tool class.
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
    /**
     *  @brief  Default constructor
     */
    ShowerStartFinderTool();

    pandora::StatusCode Run(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &peakDirection, const pandora::HitType hitType, 
        const pandora::CaloHitList &showerSpineHitList, pandora::CartesianVector &showerStartPosition, pandora::CartesianVector &showerStartDirection);

private:
    typedef std::map<const pandora::CaloHit*, float> LongitudinalPositionMap;
    typedef std::map<int, float> EnergySpectrumMap;
    typedef std::map<int, pandora::CaloHitList> LayerToHitMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create the [shower spine hit -> shower spine fit longitudinal projection] map
     *
     *  @param  spineTwoDSlidingFit the shower spine fit
     *  @param  showerSpineHitList the shower spine hit list
     *  @param  longitudinalPositionMap the output shower spine longitudinal projection map
     */
    void ObtainLongitudinalDecomposition(const TwoDSlidingFitResult &spineTwoDSlidingFit, const pandora::CaloHitList &showerSpineHitList, 
        LongitudinalPositionMap &longitudinalPositionMap) const;

    /**
     *  @brief  Create the longituidnal energy distribution 
     *
     *  @param  showerSpineHitList the shower spine hit list
     *  @param  longitudinalPositionMap the shower spine longitudinal projection map 
     *  @param  energySpectrumMap the output [longitudial projection bin -> contained energy] map
     */
    void GetEnergyDistribution(const pandora::CaloHitList &showerSpineHitList, const LongitudinalPositionMap &longitudinalPositionMap, 
        EnergySpectrumMap &energySpectrumMap) const;

    /**
     *  @brief  Find the position at which the shower cascade looks to originate, and its initial direction
     *
     *  @param  pShowerPfo the shower pfo
     *  @param  hitType the 2D view
     *  @param  spineTwoDSlidingFit the shower spine fit
     *  @param  energySpectrumMap the [longitudial projection bin -> contained energy] map 
     *  @param  showerSpineHitList the shower spine hit list
     *  @param  isEndDownstream whether the shower direction is downstream (in Z) of the neutrino vertex
     *  @param  showerStartPosition the position at which the shower cascade looks to originate
     *  @param  showerStartDirection the initial direction of the shower cascade
     */
    void FindShowerStartAndDirection(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const TwoDSlidingFitResult &spineTwoDSlidingFit, const EnergySpectrumMap &energySpectrumMap, const pandora::CaloHitList &showerSpineHitList, 
        const bool isEndDownstream, pandora::CartesianVector &showerStartPosition, pandora::CartesianVector &showerStartDirection) const;

    /**
     *  @brief  Find the mean and standard deviation of the energy depositions in the initial region
     *
     *  @param  energySpectrumMap the [longitudial projection bin -> contained energy] map
     *  @param  isEndDownstream whether the shower direction is downstream (in Z) of the neutrino vertex 
     *  @param  meanEnergy the output mean energy
     *  @param  energySigma the output standard deviation
     */
    void CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, 
        float &energySigma) const;

    /**
     *  @brief  Whether a sensible shower cascade looks to originate at a given position
     *
     *  @param  pShowerPfo the shower pfo
     *  @param  hitType the 2D view
     *  @param  spineTwoDSlidingFit the shower spine fit 
     *  @param  longitudinalDistance the longitudinal projection of the candidate shower start position
     *  @param  showerSpineHitList the shower spine hit list     
     *  @param  isEndDownstream whether the shower direction is downstream (in Z) of the neutrino vertex 
     *
     *  @return whether a sensible shower cascade looks to originate at the given position
     */
    bool IsShowerTopology(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, const TwoDSlidingFitResult &spineTwoDSlidingFit, 
        const float longitudinalDistance, const pandora::CaloHitList &showerSpineHitList, const bool isEndDownstream) const;

    /**
     *  @brief  Determine the (X,Y,Z) position and direction at a given longitudinal distance along the spine
     *
     *  @param  spineTwoDSlidingFit the shower spine fit   
     *  @param  longitudinalDistance the input longitudinal distance
     *  @param  globalPosition the output (X,Y,Z) position
     *  @param  globalDirection the output (X,Y,Z) direction
     */
    void ConvertLongitudinalProjectionToGlobal(const TwoDSlidingFitResult &spineTwoDSlidingFit, const float longitudinalDistance, 
        pandora::CartesianVector &globalPosition, pandora::CartesianVector &globalDirection) const;

    /**
     *  @brief  Build the downstream 'shower region' at a given longitudinal distance along the spine
     *
     *  @param  pShowerPfo the shower pfo 
     *  @param  hitType the 2D view
     *  @param  showerSpineHitList the shower spine hit list    
     *  @param  showerStartPosition the candidate shower start position
     *  @param  showerStartDirection the candidate shower start direction
     *  @param  isEndDownstream whether the shower direction is downstream (in Z) of the neutrino vertex
     *  @param  showerRegionPositionVector the ouput vector of shower region hit positions
     *
     *  @return whether the 'shower region finder' mechanics could proceed
     */
    pandora::StatusCode BuildShowerRegion(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, const pandora::CaloHitList &showerSpineHitList, 
        const pandora::CartesianVector &showerStartPosition, const pandora::CartesianVector &showerStartDirection, const bool isEndDownstream, 
        pandora::CartesianPointVector &showerRegionPositionVector) const;

    /**
     *  @brief  Parameterise the topological structure of the shower region
     *
     *  @param  showerRegionPositionVector the vector of shower region hit positions
     *  @param  showerStartPosition the shower start position
     *  @param  isEndDownstream whether the shower direction is downstream (in Z) of the neutrino vertex 
     *  @param  showerStartDirection the shower start direction
     *  @param  positiveEdgeStart the start position of one shower boundary   
     *  @param  positiveEdgeEnd the end position of one shower boundary
     *  @param  negativeEdgeStart the start position of the other shower boundary
     *  @param  negativeEdgeEnd the end position of the other shower boundary
     *  @param  isBetween if the shower core is between either the start or end boundary positions
     *
     *  @return whether the 'shower characterisation' mechanics could proceed
     */
    pandora::StatusCode CharacteriseShowerTopology(const pandora::CartesianPointVector &showerRegionPositionVector, const pandora::CartesianVector &showerStartPosition,
        const bool isEndDownstream, const pandora::CartesianVector &showerStartDirection, pandora::CartesianVector &positiveEdgeStart, pandora::CartesianVector &positiveEdgeEnd, 
        pandora::CartesianVector &negativeEdgeStart, pandora::CartesianVector &negativeEdgeEnd, bool &isBetween) const;

    /**
     *  @brief  Determine whether a point lies on the RHS or LHS (wrt +ve Z) of the shower core
     *
     *  @param  showerStartDirection the shower start direction
     *  @param  displacementVector the input position wrt the shower start position
     *
     *  @return whether a point lies on the RHS or LHS (wrt +ve Z) of the shower core 
     */
    bool IsClockwiseRotation(const pandora::CartesianVector &showerStartDirection, const pandora::CartesianVector &displacementVector) const;

    /**
     *  @brief  Determine the start and end positions of a shower boundary
     *
     *  @param  spineTwoDSlidingFit the shower spine fit    
     *  @param  layerFitResultMap the layer fit result map of the shower boundary fit
     *  @param  showerStartPosition the shower start position       
     *  @param  showerStartLayer the shower start layer wrt the shower region fit
     *  @param  showerEndLayer the shower end layer wrt the shower region fit
     *  @param  boundaryEdgeStart the output boundary start position
     *  @param  boundaryEdgeEnd the output boundary end position
     *
     *  @return whether suitable positions could be found
     */
    pandora::StatusCode GetBoundaryExtremalPoints(const TwoDSlidingShowerFitResult &showerTwoDSlidingFit, const LayerFitResultMap &layerFitResultMap, 
        const pandora::CartesianVector &showerStartPosition, const int showerStartLayer, const int showerEndLayer, pandora::CartesianVector &boundaryEdgeStart, 
        pandora::CartesianVector &boundaryEdgeEnd) const;

    unsigned int m_spineSlidingFitWindow;   ///< The sliding window used to fit the shower spine
    float m_longitudinalCoordinateBinSize;  ///< The longitudinal coordinate bin size 
    unsigned int m_nInitialEnergyBins;      ///< The number of longitudinal bins that define the initial region
    float m_minSigmaDeviation;              ///< The min. average energy deviation of a candidate shower start
    float m_maxEdgeGap;                     ///< The max. allowed layer gap in a shower boundary 
    float m_longitudinalShowerFraction;     ///< The shower region fraction considered
    float m_minShowerOpeningAngle;          ///< The min. opening angle of a sensible shower
    float m_molliereRadius;                 ///< The max. distance from the shower core of a collected shower region hit
    unsigned int m_showerSlidingFitWindow;  ///< The sliding window used to fit the shower region
    unsigned int m_maxLayerSeparation;      ///< The max. allowed separation between the shower start and boundary start layers
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_FINDER_TOOL_H
