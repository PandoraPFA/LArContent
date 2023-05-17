/**
 *  @file   larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h
 *
 *  @brief  Header file for the connection pathway helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CONNECTION_PATHWAY_HELPER_H
#define LAR_CONNECTION_PATHWAY_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "Objects/CartesianVector.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace pandora
{
    class CartesianVector;
    class Pandora;
} // namespace pandora

namespace lar_content
{

/**
 *  @brief  LArConnectionPathwayHelper class
 */
class LArConnectionPathwayHelper
{
public:

    enum Consistency
    {
        POSITION,
        DIRECTION,
        X_PROJECTION
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

    /**
     *  @brief  ElectronTreeVariables class
     */
    class ElectronTreeVariables
    {
    public:
        /**
         *  @brief  Constructor
         */
        ElectronTreeVariables();

        /**
         *  @brief  Assignment operator
         */
        ElectronTreeVariables &operator=(const ElectronTreeVariables &rhs);

        float m_nConnectionPathways;
        float m_showerStartX;
        float m_showerStartY;
        float m_showerStartZ;
        float m_pathwayLengthMin;
        float m_pathwayLengthMiddle;
        float m_pathwayLengthMax;
        float m_pathwayShowerStartDelta;
        float m_pathwayMaxScatteringAngle;
        float m_minShowerStartPathwayScatteringAngle2D;
        float m_middleShowerStartPathwayScatteringAngle2D;
        float m_maxShowerStartPathwayScatteringAngle2D;
        float m_pathwayEnergyMeanU;
        float m_pathwayEnergyMeanV;
        float m_pathwayEnergyMeanW;
        float m_pathwayEnergySigmaU;
        float m_pathwayEnergySigmaV;
        float m_pathwayEnergySigmaW;
        float m_pathwayMinEnergyMean;
        float m_pathwayMiddleEnergyMean;
        float m_pathwayMaxEnergyMean;
        float m_pathwayMinEnergySigma;
        float m_pathwayMiddleEnergySigma;
        float m_pathwayMaxEnergySigma;
        float m_postShowerStartLength;
        float m_postShowerStartNHits;
        float m_postShowerStartNHitsW;
        float m_minNPostShowerStartHits;
        float m_middleNPostShowerStartHits;
        float m_maxNPostShowerStartHits;
        float m_postShowerStartScatterAngle;
        float m_postShowerStartScatterAngleW;
        float m_minPostShowerStartScatterAngle;
        float m_middlePostShowerStartScatterAngle;
        float m_maxPostShowerStartScatterAngle;
        float m_postShowerStartOpeningAngleW;
        float m_minPostShowerStartOpeningAngle;
        float m_middlePostShowerStartOpeningAngle;
        float m_maxPostShowerStartOpeningAngle;
        float m_postShowerStartOpeningAngleAsymmetryW;
        float m_minPostShowerStartOpeningAngleAsymmetry;
        float m_middlePostShowerStartOpeningAngleAsymmetry;
        float m_maxPostShowerStartOpeningAngleAsymmetry;
        float m_postShowerStartNuVertexHitAsymmetryW;
        float m_minPostShowerStartNuVertexHitAsymmetry;
        float m_middlePostShowerStartNuVertexHitAsymmetry;
        float m_maxPostShowerStartNuVertexHitAsymmetry;
        float m_postShowerStartNuVertexEnergyAsymmetryW;
        float m_minPostShowerStartNuVertexEnergyAsymmetry;
        float m_middlePostShowerStartNuVertexEnergyAsymmetry;
        float m_maxPostShowerStartNuVertexEnergyAsymmetry;
        float m_postShowerStartShowerStartHitAsymmetryW;
        float m_minPostShowerStartShowerStartHitAsymmetry;
        float m_middlePostShowerStartShowerStartHitAsymmetry;
        float m_maxPostShowerStartShowerStartHitAsymmetry;
        float m_postShowerStartShowerStartEnergyAsymmetryW;
        float m_minPostShowerStartShowerStartEnergyAsymmetry;
        float m_middlePostShowerStartShowerStartEnergyAsymmetry;
        float m_maxPostShowerStartShowerStartEnergyAsymmetry;        
        float m_postShowerStartNuVertexMeanRadialDistanceW;
        float m_minPostShowerStartNuVertexMeanRadialDistance;
        float m_middlePostShowerStartNuVertexMeanRadialDistance;
        float m_maxPostShowerStartNuVertexMeanRadialDistance;
        float m_postShowerStartNuVertexEnergyWeightedMeanRadialDistanceW;
        float m_minPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
        float m_middlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
        float m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
        float m_postShowerStartShowerStartMeanRadialDistanceW;
        float m_minPostShowerStartShowerStartMeanRadialDistance;
        float m_middlePostShowerStartShowerStartMeanRadialDistance;
        float m_maxPostShowerStartShowerStartMeanRadialDistance;
        float m_postShowerStartShowerStartEnergyWeightedMeanRadialDistanceW;
        float m_minPostShowerStartShowerStartEnergyWeightedMeanRadialDistance;
        float m_middlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance;
        float m_maxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance;
        float m_postShowerStartNuVertexMoliereRadiusW;
        float m_minPostShowerStartNuVertexMoliereRadius;
        float m_middlePostShowerStartNuVertexMoliereRadius;
        float m_maxPostShowerStartNuVertexMoliereRadius;
        float m_postShowerStartShowerStartMoliereRadiusW;
        float m_minPostShowerStartShowerStartMoliereRadius;
        float m_middlePostShowerStartShowerStartMoliereRadius;
        float m_maxPostShowerStartShowerStartMoliereRadius;
        float m_positiveOpeningAngleW;
        float m_negativeOpeningAngleW;
        float m_maxOpeningAngleW;
        float m_showerApexLW;
        float m_minShowerApexL;
        float m_middleShowerApexL;
        float m_maxShowerApexL;
        float m_showerApexTW;
        float m_minShowerApexT;
        float m_middleShowerApexT;
        float m_maxShowerApexT;
        float m_foundHitRatioW;
        float m_minFoundHitRatio;
        float m_middleFoundHitRatio;
        float m_maxFoundHitRatio;
        float m_fitShowerStartLW;
        float m_fitShowerStartTW;
        float m_postShowerStartMinHalfOpeningAngle;
        float m_postShowerStartMaxHalfOpeningAngle;
        float m_postShowerStartOpeningAngle;
        float m_postShowerStartOpeningAngleAsymmetry;
        float m_postShowerStartMeanTransverseAngle;
        float m_postShowerStartMeanLWeightedTransverseAngle;
        float m_postShowerStartMeanRadialDistance;
        float m_postShowerStartRadialDistanceSigma;
        float m_postShowerStartEnergyWeightedMeanRadialDistance;
        float m_postShowerStartEstimatedMoliereRadius;
        float m_postShowerStartLWeightedMeanRadialDistance;
        float m_postShowerStartLWeightedRadialDistanceSigma;
        float m_postShowerStartInitialGapSize;
        float m_postShowerStartMaxGapSize;
        float m_initialRegionDistanceToNuVertex;
        float m_initialRegionDistanceInGaps;
        float m_initialRegionMaxGapSize;
        float m_initialGapSizeW;
        float m_minInitialGapSize;
        float m_middleInitialGapSize;
        float m_maxInitialGapSize;
        float m_largestGapSizeW;
        float m_minLargestGapSize;
        float m_middleLargestGapSize;
        float m_maxLargestGapSize;
        float m_largestProjectedGapSizeW;
        float m_minLargestProjectedGapSize;
        float m_middleLargestProjectedGapSize;
        float m_maxLargestProjectedGapSize;
        float m_hitLineDensityW;
        float m_minHitLineDensity;
        float m_middleHitLineDensity;
        float m_maxHitLineDensity;
        float m_nViewsWithAmbiguousHits;
        float m_nAmbiguousHits2D;
        float m_minNAmbiguousHits;
        float m_maxNAmbiguousHits;
        float m_ambiguousHitUnaccountedEnergyU;
        float m_ambiguousHitUnaccountedEnergyV;
        float m_ambiguousHitUnaccountedEnergyW;
        float m_ambiguousHitMinUnaccountedEnergy;
        float m_ambiguousHitMaxUnaccountedEnergy;
        float m_ambiguousHitShowerEnergyRatioU;
        float m_ambiguousHitShowerEnergyRatioV;
        float m_ambiguousHitShowerEnergyRatioW;
        float m_ambiguousHitMinShowerEnergyRatio;
        float m_ambiguousHitMaxShowerEnergyRatio;
    };

    static void FillElectronTreeVariables(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const pandora::CartesianVector &nuVertexPosition, 
        const pandora::CaloHitList *const pCaloHitListU, const pandora::CaloHitList *const pCaloHitListV, const pandora::CaloHitList *const pCaloHitListW,
        const LArConnectionPathwayHelper::Consistency &consistency, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool FillAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const pandora::CartesianVector &nuVertexPosition,
        pandora::CartesianVector &middleShowerStart3D, const pandora::CaloHitList *const pCaloHitListU, const pandora::CaloHitList *const pCaloHitListV, 
        const pandora::CaloHitList *const pCaloHitListW, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool GetViewAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
        const pandora::HitType hitType, const pandora::CartesianVector &nuVertexPosition, const pandora::CartesianVector &connectionPathwayDirection3D,
        const pandora::CaloHitList *const pCaloHitList, float &unaccountedHitEnergy, float &showerEnergyRatio);

    static pandora::CaloHitList FindAmbiguousContinuousSpine(const pandora::CaloHitList &caloHitList, const pandora::CaloHitList &ambiguousHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition);

    static bool FillInitialRegionVariables(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, 
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW,
        pandora::CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool FindPostShowerStart2DVariables(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, const ProtoShower &protoShower, 
        const pandora::HitType hitType, float &initialGapSize, float &maxGapSize, float &maxProjectedGapSize, float &hitLineDensity);

    static bool FillPostShowerStartVariables(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, 
        const ProtoShower &protoShowerW, pandora::CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool FindPostShowerStart2DVariables(pandora::Algorithm *const pAlgorithm, const pandora::CaloHitList &spineHitList, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::HitType hitType, const bool isDownstream, const pandora::CartesianVector &projectedNuVertex, const pandora::CartesianVector &projectedShowerStart, 
        const pandora::CartesianVector &initialShowerDirection, const pandora::CartesianVector &connectionPathwayDirection, int &nHits, float &scatterAngle,
        float &openingAngle, float &openingAngleAsymmetry, float &nuVertexHitAsymmetry, float &nuVertexEnergyAsymmetry, float &showerStartHitAsymmetry, 
        float &showerStartEnergyAsymmetry, float &nuVertexMeanRadialDistance, float &nuVertexEnergyWeightedMeanRadialDistance, float &showerStartMeanRadialDistance, 
        float &showerStartEnergyWeightedRadialDistance, float &nuVertexMoliereRadius, float &showerStartMoliereRadius, float &positiveOpeningAngle, 
        float &negativeOpeningAngle, float &showerApexL, float &showerApexT, float &fitShowerStartL, float &fitShowerStartT, float &foundHitRatio);

    static bool FindShowerStarts3D(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const LArConnectionPathwayHelper::Consistency &consistency, 
        const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &minShowerStart3D, pandora::CartesianVector &middleShowerStart3D, 
        pandora::CartesianVector &maxShowerStart3D);

    static bool FillConnectionPathwayVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &minShowerStart3D,
        pandora::CartesianVector &middleShowerStart3D, pandora::CartesianVector &maxShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static float GetLargest3DKink(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFitU, const TwoDSlidingFitResult &spineFitV, 
        const TwoDSlidingFitResult &spineFitW, const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &maxShowerStart3D);

    static float GetLargest3DKinkFromView(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, 
        const TwoDSlidingFitResult &spineFit1, const TwoDSlidingFitResult &spineFit2, const pandora::HitType hitType, const pandora::HitType hitType1, 
        const pandora::HitType hitType2, const pandora::CartesianVector &maxShowerStart3D);

    static float Get2DKink(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFitU,
        const TwoDSlidingFitResult &spineFitV, const TwoDSlidingFitResult &spineFitW, const pandora::CartesianVector &showerStart3D);

    static float GetLargest2DKinkFromView(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, 
        const pandora::HitType hitType, const pandora::CartesianVector &showerStart3D);

    static bool FillConnectionPathwayEnergyVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &middleShowerStart3D, 
        LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static pandora::CaloHitList GetConnectionPathwayProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
        const pandora::CartesianVector &showerStart3D, const pandora::HitType hitType);

    static float CharacteriseEnergyMean(const pandora::CaloHitList &caloHitList);

    static float CharacteriseEnergySigma(const pandora::CaloHitList &caloHitList, const float energyMean);


   static bool AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation);

   static bool AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector &showerStartU, 
       const pandora::CartesianVector &showerStartV, const pandora::CartesianVector &showerStartW, const float maxXSeparation, const float maxSeparation, float &metric);

   static bool AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
      const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle, float &metric);

   static bool AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle);

   static bool AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, pandora::CartesianVector directionU,
       pandora::CartesianVector directionV, pandora::CartesianVector directionW, const float maxOpeningAngle, float &metric);

   static bool FindDirection3D(pandora::Algorithm *const pAlgorithm, const bool isDownstream, const pandora::CartesianVector &startPositionU, 
       const pandora::CartesianVector &startPositionV, const pandora::CartesianVector &startPositionW, const pandora::CartesianVector &startPosition3D, 
       const pandora::CartesianVector &directionU, const pandora::CartesianVector &directionV, const pandora::CartesianVector &directionW, pandora::CartesianVector &direction3D);

   static bool FindShowerVertexFromPosition(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

    static bool FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &uShowerStart3D, 
        pandora::CartesianVector &vShowerStart3D, pandora::CartesianVector &wShowerStart3D);

   static bool FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector nuVertexPosition, 
       const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

   static bool FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector nuVertexPosition, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
        const float maxSeparation, pandora::CartesianVector &showerStart3D);

    static bool FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, 
        pandora::CartesianVector &showerStart3D);

    static bool FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ProtoShower &protoShower, const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, pandora::CartesianVector &showerStart3D);

    static void GetMinMiddleMax(const float value1, const float value2, const float value3, float &minValue, float &middleValue,
        float &maxValue);
};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_HELPER_H
