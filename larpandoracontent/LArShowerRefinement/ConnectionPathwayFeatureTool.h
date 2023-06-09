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

    typedef MvaFeatureTool<const pandora::Algorithm *const, const pandora::ParticleFlowObject *const, const pandora::CartesianVector&, const ProtoShowerMatch&, 
        const pandora::CartesianPointVector&> ConnectionPathwayFeatureTool;

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

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetViewInitialRegionVariables(const pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, 
        const pandora::HitType hitType, float &initialGapSize, float &largestGapSize);

    unsigned int m_nHitsToConsider;
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

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float Get2DKink(const pandora::Algorithm *const pAlgorithm, const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianVector &showerStart3D) const;
    float GetLargest2DKinkFromView(const pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, const pandora::HitType hitType, const pandora::CartesianVector &showerStart3D) const;

    unsigned int m_spineFitWindow;
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

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetViewShowerRegionVariables(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const ProtoShowerMatch &protoShowerMatch, const pandora::HitType hitType, const pandora::CartesianVector &showerStart3D, float &nHits, float &foundHitRatio, float &scatterAngle, 
        float &openingAngle, float &nuVertexEnergyAsymmetry, float &nuVertexEnergyWeightedMeanRadialDistance, float &showerStartEnergyAsymmetry, float &showerStartMoliereRadius);

    void BuildViewShower(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const TwoDSlidingFitResult &spineFit,
        const pandora::HitType hitType, const pandora::CartesianVector &showerStart2D, const pandora::CartesianVector &nuVertex2D, pandora::CaloHitList &postShowerHitList, 
        pandora::CartesianPointVector &postShowerPositions);

    void GetShowerHitVariables(const pandora::CaloHitList &spineHitList, const pandora::CaloHitList &postShowerHitList, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::HitType hitType, float &nHits, float &foundHitRatio);

    void CalculateViewScatterAngle(const pandora::CartesianVector &nuVertex2D, const TwoDSlidingFitResult &spineFitResult, const pandora::CartesianVector &showerStart2D, 
        const TwoDSlidingFitResult &showerFitResult, float &scatterAngle);

    void CalculateViewOpeningAngle(const pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &showerFitResult, 
        const pandora::CaloHitList &postShowerHitList, const pandora::CartesianVector &showerStart2D, float &openingAngle);

    void CalculateViewNuVertexConsistencyVariables(const TwoDSlidingFitResult &spineFitResult, const pandora::CaloHitList &postShowerHitList, 
        const bool isDownstream, const pandora::CartesianVector &nuVertex2D, float &nuVertexEnergyAsymmetry, float &nuVertexEnergyWeightedMeanRadialDistance);

    void CalculateViewShowerStartConsistencyVariables(const TwoDSlidingFitResult &showerFitResult, const pandora::CaloHitList &postShowerHitList, 
        const bool isShowerDownstream, float &showerStartEnergyAsymmetry, float &showerStartMoliereRadius);

    unsigned int m_spineFitWindow;
    float m_showerRadius;
    unsigned int m_showerFitWindow;
    float m_edgeStep;
    float m_moliereFraction;
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

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianPointVector &showerStarts3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void CalculateNAmbiguousViews(const ProtoShowerMatch &protoShowerMatch, float &nAmbiguousViews);

    bool GetViewAmbiguousHitVariables(const pandora::Algorithm *const pAlgorithm, const ProtoShowerMatch &protoShowerMatch, 
        const pandora::HitType hitType, const pandora::CartesianVector &nuVertex3D, float &unaccountedHitEnergy);

    void BuildAmbiguousSpines(const pandora::Algorithm *const pAlgorithm, const pandora::HitType hitType, const ProtoShower &protoShower, 
        const pandora::CartesianVector &nuVertex2D, std::map<int, pandora::CaloHitList> &ambiguousHitSpines, pandora::CaloHitList &hitsToExcludeInEnergyCalcs);

    pandora::StatusCode GetHitListOfType(const pandora::Algorithm *const pAlgorithm, const pandora::HitType hitType, const pandora::CaloHitList *&pCaloHitList) const;

    pandora::CaloHitList FindAmbiguousContinuousSpine(const pandora::CaloHitList &caloHitList, const pandora::CaloHitList &ambiguousHitList, 
        const pandora::CartesianVector &nuVertex2D);

    std::string m_caloHitListNameU;
    std::string m_caloHitListNameV;
    std::string m_caloHitListNameW;
    float m_maxTransverseDistance;
    unsigned int m_maxSampleHits;
    float m_maxHitSeparation;
    float m_maxTrackFraction;
};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_FEATURE_TOOLS_H
