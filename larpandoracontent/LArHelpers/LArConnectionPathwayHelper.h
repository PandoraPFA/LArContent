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

        float m_pathwayLength;
        float m_pathwayShowerStartDelta;
        float m_pathwayMaxScatteringAngle;
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
        int m_postShowerStartNHits;
        float m_postShowerStartScatterAngle;
        float m_postShowerStartMeanTransverseAngle;
        float m_postShowerStartMeanRadialDistance;
        float m_postShowerStartRadialDistanceSigma;
        float m_postShowerStartEnergyWeightedMeanRadialDistance;
        float m_postShowerStartEstimatedMoliereRadius;
        float m_postShowerStartInitialGapSize;
        float m_postShowerStartMaxGapSize;
        float m_initialRegionDistanceToNuVertex;
        float m_initialRegionDistanceInGaps;
        float m_initialRegionMaxGapSize;
    };

    static void FillElectronTreeVariables(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const pandora::CartesianVector &nuVertexPosition, 
        const LArConnectionPathwayHelper::Consistency &consistency, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool FillAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW, const pandora::CartesianVector &nuVertexPosition,
        pandora::CartesianVector &middleShowerStart3D, const TwoDSlidingFitResult &spineFitU, const TwoDSlidingFitResult &spineFitV, const TwoDSlidingFitResult &spineFitW,
        const pandora::CaloHitList *const pCaloHitListU, const pandora::CaloHitList *const pCaloHitListV, const pandora::CaloHitList *const pCaloHitListW);

    static void GetViewAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const ElectronProtoShower &protoShower, 
        const pandora::HitType hitType, const pandora::CartesianVector &nuVertexPosition, const pandora::CartesianVector &connectionPathwayDirection3D,
        const pandora::CaloHitList *const pCaloHitList, float &unaccountedHitEnergy, float &showerEnergyRatio);

    static pandora::CaloHitList FindAmbiguousContinuousSpine(const pandora::CaloHitList &caloHitList, const pandora::CaloHitList &ambiguousHitList, 
        const pandora::CartesianVector &projectedNuVertexPosition);

    static bool FillInitialRegionVariables(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, 
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW,
        pandora::CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool FillPostShowerStartVariables(pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const pandora::CartesianVector &nuVertexPosition, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, 
        const ProtoShower &protoShowerW, pandora::CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static bool FindShowerStarts3D(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const LArConnectionPathwayHelper::Consistency &consistency, 
        const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &minShowerStart3D, pandora::CartesianVector &middleShowerStart3D, 
        pandora::CartesianVector &maxShowerStart3D);

    static bool FillConnectionPathwayVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &minShowerStart3D,
        pandora::CartesianVector &middleShowerStart3D, pandora::CartesianVector &maxShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    static float GetLargest3DKink(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, 
        const ProtoShower &protoShowerW, const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &maxShowerStart3D);

    static float GetLargest3DKinkFromView(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, 
        const TwoDSlidingFitResult &spineFit1, const TwoDSlidingFitResult &spineFit2, const pandora::HitType hitType, const pandora::HitType hitType1, 
        const pandora::HitType hitType2, const pandora::CartesianVector &maxShowerStart3D);

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

    static bool FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector nuVertexPosition,
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, 
        pandora::CartesianVector &showerStart3D);

    static bool FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, 
        pandora::CartesianVector &showerStart3D);

};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_HELPER_H
