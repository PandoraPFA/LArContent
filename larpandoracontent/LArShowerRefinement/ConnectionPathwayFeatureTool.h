/**
 *  @file   larpandoracontent/LArShowerRefinement/ConnectionPathwayFeatureTool.h
 *
 *  @brief  Header file for the connection pathway feature tools
 *
 *  $Log: $
 */
#ifndef LAR_CONNECTION_PATHWAY_FEATURE_TOOLS_H
#define LAR_CONNECTION_PATHWAY_FEATURE_TOOLS_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

namespace lar_content
{

typedef MvaFeatureTool<const pandora::Algorithm *const, const pandora::ParticleFlowObject *const, const pandora::CartesianVector &, const ProtoShowerMatch &, const pandora::CartesianPointVector &> ConnectionPathwayFeatureTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  InitialRegionFeatureTool to calculate variables related to the initial shower region
 */
class InitialRegionFeatureTool : public ConnectionPathwayFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    InitialRegionFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculate the initial region variables for the input view
     *
     *  @param  pAlgorithm the algorithm
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  protoShowerMatch the ProtoShower match
     *  @param  hitType the 2D view
     *  @param  initialGapSize the output initial gap size
     *  @param  largestGapSize the output largest gap size
     */
    void GetViewInitialRegionVariables(const pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::HitType hitType, float &initialGapSize, float &largestGapSize) const;

    float m_defaultFloat;           ///< Default float value
    unsigned int m_nHitsToConsider; ///< The number of hits which defines the initial region
    float m_maxInitialGapSizeLimit; ///< maxInitialGapSizeLimit max. limit
    float m_minLargestGapSizeLimit; ///< minLargestGapSizeLimit max. limit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ConnectionRegionFeatureTool to calculate variables related to the connection pathway region
 */
class ConnectionRegionFeatureTool : public ConnectionPathwayFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ConnectionRegionFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Obtain a cautious estimate of the largest 2D deflection of the connection pathway
     *
     *  @param  pAlgorithm the algorithm
     *  @param  protoShowerMatch the ProtoShower match
     *  @param  showerStart3D the 3D shower start position
     *
     *  @return  a cautious estimate of the largest 2D deflection of the connection pathway 
     */
    float Get2DKink(const pandora::Algorithm *const pAlgorithm, const ProtoShowerMatch &protoShowerMatch,
        const pandora::CartesianVector &showerStart3D) const;

    /**
     *  @brief  Obtain a cautious estimate of the largest 2D deflection of a connection pathway in a given view
     *
     *  @param  pAlgorithm the algorithm
     *  @param  spineFit the shower spine fit
     *  @param  hitType the 2D view  
     *  @param  showerStart3D the 3D shower start position
     *
     *  @return a cautious estimate of the largest 2D deflection of a connection pathway in a given view
     */
    float GetLargest2DKinkFromView(const pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit,
        const pandora::HitType hitType, const pandora::CartesianVector &showerStart3D) const;

    float m_defaultFloat;                  ///< Default float value
    unsigned int m_spineFitWindow;         ///< The spine fit window
    int m_nLayersHalfWindow;               ///< The half window of each segment
    float m_pathwayLengthLimit;            ///< pathwayLengthLimit max. limit
    float m_pathwayScatteringAngle2DLimit; ///< pathwayScatteringAngle2DLimit max. limit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ShowerRegionFeatureTool to calculate variables related to the shower region
 */
class ShowerRegionFeatureTool : public ConnectionPathwayFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerRegionFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculate the shower region variables for the input view
     *
     *  @param  pAlgorithm the algorithm
     *  @param  pShowerPfo the shower pfo
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  protoShowerMatch the ProtoShower match
     *  @param  hitType the 2D view
     *  @param  showerStart3D the 3D shower start position
     *  @param  nHits the output number of shower region hits
     *  @param  foundHitRatio the output found hit ratio
     *  @param  scatterAngle the output scatter angle
     *  @param  openingAngle the output opening angle
     *  @param  nuVertexEnergyAsymmetry the output neutrino vertex energy asymmetry
     *  @param  nuVertexEnergyWeightedMeanRadialDistance the output neutrino vertex energy weighted mean radial distance
     *  @param  showerStartEnergyAsymmetry the output shower start energy asymmetry 
     *  @param  showerStartMoliereRadius the output shower start moliere radius
     */
    void GetViewShowerRegionVariables(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const pandora::HitType hitType,
        const pandora::CartesianVector &showerStart3D, float &nHits, float &foundHitRatio, float &scatterAngle, float &openingAngle,
        float &nuVertexEnergyAsymmetry, float &nuVertexEnergyWeightedMeanRadialDistance, float &showerStartEnergyAsymmetry,
        float &showerStartMoliereRadius);

    /**
     *  @brief  Collect the shower region hits in a given view
     *
     *  @param  pShowerPfo the shower pfo
     *  @param  spineFit the shower spine fit
     *  @param  hitType the 2D view
     *  @param  showerStart2D the 2D shower start position
     *  @param  nuVertex2D the 2D neutrino vertex position
     *  @param  postShowerHitList the collected shower region hit list
     *  @param  postShowerPositions the collected shower region hit position vector
     */
    void BuildViewShower(const pandora::ParticleFlowObject *const pShowerPfo, const TwoDSlidingFitResult &spineFit,
        const pandora::HitType hitType, const pandora::CartesianVector &showerStart2D, const pandora::CartesianVector &nuVertex2D,
        pandora::CaloHitList &postShowerHitList, pandora::CartesianPointVector &postShowerPositions);

    /**
     *  @brief  Evaluate the variables associated with the shower region hit multiplicity
     *
     *  @param  spineHitList the shower spine hit list
     *  @param  postShowerHitList the collected shower region hit list
     *  @param  pShowerPfo the shower pfo
     *  @param  hitType the 2D view
     *  @param  nHits the output number of shower region hits
     *  @param  foundHitRatio the output found hit ratio
     */
    void GetShowerHitVariables(const pandora::CaloHitList &spineHitList, const pandora::CaloHitList &postShowerHitList,
        const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, float &nHits, float &foundHitRatio);

    /**
     *  @brief  Calculate the connection pathway-shower region scatter angle
     *
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  spineFitResult the shower spine fit
     *  @param  showerStart2D the 2D shower start position
     *  @param  showerFitResult the shower region fit
     *  @param  scatterAngle the output scatter angle
     */
    void CalculateViewScatterAngle(const pandora::CartesianVector &nuVertex2D, const TwoDSlidingFitResult &spineFitResult,
        const pandora::CartesianVector &showerStart2D, const TwoDSlidingFitResult &showerFitResult, float &scatterAngle);

    /**
     *  @brief  Calculate the opening angle of the shower region
     *
     *  @param  showerFitResult the shower region fit
     *  @param  postShowerHitList the collected shower region hit list 
     *  @param  showerStart2D the 2D shower start position 
     *  @param  openingAngle the output opening angle
     */
    void CalculateViewOpeningAngle(const TwoDSlidingFitResult &showerFitResult, const pandora::CaloHitList &postShowerHitList, 
        const pandora::CartesianVector &showerStart2D, float &openingAngle);

    /**
     *  @brief  Evaluate the neutrino vertex consistency variables
     *
     *  @param  spineHitList the shower spine hit list
     *  @param  postShowerHitList the collected shower region hit list
     *  @param  isDownstream whether the shower direction is downstream (in Z) of the neutrino vertex
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  nuVertexEnergyAsymmetry the output neutrino vertex energy asymmetry
     *  @param  nuVertexEnergyWeightedMeanRadialDistance the output neutrino vertex energy weighted mean radial distance
     */
    void CalculateViewNuVertexConsistencyVariables(const TwoDSlidingFitResult &spineFitResult,
        const pandora::CaloHitList &postShowerHitList, const bool isDownstream, const pandora::CartesianVector &nuVertex2D,
        float &nuVertexEnergyAsymmetry, float &nuVertexEnergyWeightedMeanRadialDistance);

    /**
     *  @brief  Evaluate the shower start consistency variables
     *
     *  @param  showerFitResult the shower fit result
     *  @param  postShowerHitList the collected shower region hit list
     *  @param  isShowerDownstream whether the shower direction is downstream (in Z) of the neutrino vertex
     *  @param  showerStartEnergyAsymmetry the output shower start energy asymmetry 
     *  @param  showerStartMoliereRadius the output shower start moliere radius
     */
    void CalculateViewShowerStartConsistencyVariables(const TwoDSlidingFitResult &showerFitResult, const pandora::CaloHitList &postShowerHitList,
        const bool isShowerDownstream, float &showerStartEnergyAsymmetry, float &showerStartMoliereRadius);

    float m_defaultFloat;           ///< Default float value
    float m_defaultRatio;           ///< Default float value for ratios
    unsigned int m_spineFitWindow;  ///< The spine fit window
    float m_showerRadius;           ///< The max. separation distance between a shower region hit and the shower core
    unsigned int m_showerFitWindow; ///< The shower fit window
    float m_edgeStep;               ///< The binning of the shower boundaries
    float m_moliereFraction;        ///< The energy fraction which corresponds to minShowerStartMoliereRadius
    float m_maxNHitsLimit;          ///< maxNHits max. limit
    float m_maxFoundHitRatioLimit;  ///< maxFoundHitRatio max. limit
    float m_maxScatterAngleLimit;   ///< maxScatterAngle max. limit
    float m_maxOpeningAngleLimit;   ///< maxOpeningAngle max. limit
    float m_maxNuVertexEnergyWeightedMeanRadialDistanceLimit; ///< maxNuVertexEnergyWeightedMeanRadialDistance max. limit
    float m_minShowerStartMoliereRadiusLimit;                 ///< minShowerStartMoliereRadius max. limit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  AmbiguousRegionFeatureTool to calculate variables related to the shower region
 */
class AmbiguousRegionFeatureTool : public ConnectionPathwayFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    AmbiguousRegionFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm,
        const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Count the number of views with ambiguous hits
     *
     *  @param  protoShowerMatch the ProtoShower match
     *  @param  nAmbiguousViews the output number of ambiguous views
     */
    void CalculateNAmbiguousViews(const ProtoShowerMatch &protoShowerMatch, float &nAmbiguousViews);

    /**
     *  @brief  Calculate the ambiguous region variables for the input view
     *
     *  @param  pAlgorithm the algorithm
     *  @param  protoShowerMatch the ProtoShower match 
     *  @param  hitType the 2D view
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  unaccountedHitEnergy the output unaccounted hit energy 
     *
     *  @return whether the ambiguous region variables could be calculated
     */
    bool GetViewAmbiguousHitVariables(const pandora::Algorithm *const pAlgorithm, const ProtoShowerMatch &protoShowerMatch,
        const pandora::HitType hitType, const pandora::CartesianVector &nuVertex3D, float &unaccountedHitEnergy);

    /**
     *  @brief  Determine the spine hits of the particles with which the ambiguous hits are shared
     *
     *  @param  pAlgorithm the algorithm
     *  @param  hitType the 2D view
     *  @param  protoShower the ProtoShower
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  ambiguousHitSpines the output [particle index -> shower spine hits] map
     *  @param  hitsToExcludeInEnergyCalcs the list of hits to exclude in energy calculations
     */
    void BuildAmbiguousSpines(const pandora::Algorithm *const pAlgorithm, const pandora::HitType hitType, const ProtoShower &protoShower,
        const pandora::CartesianVector &nuVertex2D, std::map<int, pandora::CaloHitList> &ambiguousHitSpines,
        pandora::CaloHitList &hitsToExcludeInEnergyCalcs);

    /**
     *  @brief  Obtain the event hit list of a given view
     *
     *  @param  pAlgorithm the algorithm 
     *  @param  hitType the 2D view 
     *  @param  pCaloHitList the output 2D hit list
     *
     *  @return whether a valid 2D hit list could be found
     */
    pandora::StatusCode GetHitListOfType(
        const pandora::Algorithm *const pAlgorithm, const pandora::HitType hitType, const pandora::CaloHitList *&pCaloHitList) const;

    /**
     *  @brief  Determine a continuous pathway of an ambigous particle's spine hits 
     *
     *  @param  caloHitList the input ambiguous particle spine hit list
     *  @param  ambiguousHitList the ambiguous hit list
     *  @param  nuVertex2D the 2D neutrino vertex
     *
     *  @return a continuous hit pathway
     */
    pandora::CaloHitList FindAmbiguousContinuousSpine(
        const pandora::CaloHitList &caloHitList, const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianVector &nuVertex2D);

    float m_defaultFloat;           ///< Default float value
    std::string m_caloHitListNameU; ///< The event U view hit list
    std::string m_caloHitListNameV; ///< The event V view hit list
    std::string m_caloHitListNameW; ///< The event W view hit list
    float m_maxTransverseDistance;  ///< The max. proximity of a hits, included in a trajectory energy calcs.
    unsigned int m_maxSampleHits;   ///< The max. number of hits considered in the spine energy calcs.
    float m_maxHitSeparation;       ///< The max. separation of connected hits
    float m_maxTrackFraction;       ///< The fraction of found hits which are considered in the energy calcs.
};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_FEATURE_TOOLS_H
