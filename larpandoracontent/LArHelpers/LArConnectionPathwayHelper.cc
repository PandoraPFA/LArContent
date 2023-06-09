/**<
 *  @file   larpandoracontent/LArHelpers/LArConnectionPathwayHelper.cc
 *
 *  @brief  Implementation of the connection pathway helper class.
 *
 *  $Log: $
 */

#include "Objects/CartesianVector.h"
#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"
#include "Pandora/Pandora.h"
#include "Pandora/Algorithm.h"

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"


#include <limits>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------
// Electron Tree
//------------------------------------------------------------------------------------------------------------------------------------------

LArConnectionPathwayHelper::ElectronTreeVariables::ElectronTreeVariables() :
    m_nConnectionPathways(-10.f),
    m_showerStartX(-9999.f),
    m_showerStartY(-9999.f),
    m_showerStartZ(-9999.f),
    m_pathwayLengthMin(-10.f),
    m_pathwayLengthMiddle(-10.f),
    m_pathwayLengthMax(-10.f),
    m_pathwayShowerStartDelta(-10.f),
    m_pathwayMaxScatteringAngle(-10.f),
    m_minShowerStartPathwayScatteringAngle2D(-10.f),
    m_middleShowerStartPathwayScatteringAngle2D(-10.f),
    m_maxShowerStartPathwayScatteringAngle2D(-10.f),
    m_pathwayEnergyMeanU(-10.f),
    m_pathwayEnergyMeanV(-10.f),
    m_pathwayEnergyMeanW(-10.f),
    m_pathwayEnergySigmaU(-10.f),
    m_pathwayEnergySigmaV(-10.f),
    m_pathwayEnergySigmaW(-10.f),
    m_pathwayMinEnergyMean(-10.f),
    m_pathwayMiddleEnergyMean(-10.f),
    m_pathwayMaxEnergyMean(-10.f),
    m_pathwayMinEnergySigma(-10.f),
    m_pathwayMiddleEnergySigma(-10.f),
    m_pathwayMaxEnergySigma(-10.f),
    m_postShowerStartLength(-10.f),
    m_postShowerStartNHits(-10),
    m_postShowerStartNHitsW(-10.f),
    m_minNPostShowerStartHits(-10.f),
    m_middleNPostShowerStartHits(-10.f),
    m_maxNPostShowerStartHits(-10.f),
    m_postShowerStartScatterAngle(-10.f),
    m_postShowerStartScatterAngleW(-10.f),
    m_minPostShowerStartScatterAngle(-10.f),
    m_middlePostShowerStartScatterAngle(-10.f),
    m_maxPostShowerStartScatterAngle(-10.f),
    m_postShowerStartOpeningAngleW(-10.f),
    m_minPostShowerStartOpeningAngle(-10.f),
    m_middlePostShowerStartOpeningAngle(-10.f),
    m_maxPostShowerStartOpeningAngle(-10.f),
    m_postShowerStartOpeningAngleAsymmetryW(-10.f),
    m_minPostShowerStartOpeningAngleAsymmetry(-10.f),
    m_middlePostShowerStartOpeningAngleAsymmetry(-10.f),
    m_maxPostShowerStartOpeningAngleAsymmetry(-10.f),
    m_postShowerStartNuVertexHitAsymmetryW(-10.f),
    m_minPostShowerStartNuVertexHitAsymmetry(-10.f),
    m_middlePostShowerStartNuVertexHitAsymmetry(-10.f),
    m_maxPostShowerStartNuVertexHitAsymmetry(-10.f),
    m_postShowerStartNuVertexEnergyAsymmetryW(-0.5f),
    m_minPostShowerStartNuVertexEnergyAsymmetry(-0.5f),
    m_middlePostShowerStartNuVertexEnergyAsymmetry(-0.5f),
    m_maxPostShowerStartNuVertexEnergyAsymmetry(-0.5f),
    m_postShowerStartShowerStartHitAsymmetryW(-10.f),
    m_minPostShowerStartShowerStartHitAsymmetry(-10.f),
    m_middlePostShowerStartShowerStartHitAsymmetry(-10.f),
    m_maxPostShowerStartShowerStartHitAsymmetry(-10.f),
    m_postShowerStartShowerStartEnergyAsymmetryW(-0.5f),
    m_minPostShowerStartShowerStartEnergyAsymmetry(-0.5f),
    m_middlePostShowerStartShowerStartEnergyAsymmetry(-0.5f),
    m_maxPostShowerStartShowerStartEnergyAsymmetry(-0.5f),        
    m_postShowerStartNuVertexMeanRadialDistanceW(-10.f),
    m_minPostShowerStartNuVertexMeanRadialDistance(-10.f),
    m_middlePostShowerStartNuVertexMeanRadialDistance(-10.f),
    m_maxPostShowerStartNuVertexMeanRadialDistance(-10.f),
    m_postShowerStartNuVertexEnergyWeightedMeanRadialDistanceW(-10.f),
    m_minPostShowerStartNuVertexEnergyWeightedMeanRadialDistance(-10.f),
    m_middlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance(-10.f),
    m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance(-10.f),
    m_postShowerStartShowerStartMeanRadialDistanceW(-10.f),
    m_minPostShowerStartShowerStartMeanRadialDistance(-10.f),
    m_middlePostShowerStartShowerStartMeanRadialDistance(-10.f),
    m_maxPostShowerStartShowerStartMeanRadialDistance(-10.f),
    m_postShowerStartShowerStartEnergyWeightedMeanRadialDistanceW(-10.f),
    m_minPostShowerStartShowerStartEnergyWeightedMeanRadialDistance(-10.f),
    m_middlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance(-10.f),
    m_maxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance(-10.f),
    m_postShowerStartNuVertexMoliereRadiusW(-10.f),
    m_minPostShowerStartNuVertexMoliereRadius(-10.f),
    m_middlePostShowerStartNuVertexMoliereRadius(-10.f),
    m_maxPostShowerStartNuVertexMoliereRadius(-10.f),
    m_postShowerStartShowerStartMoliereRadiusW(-10.f),
    m_minPostShowerStartShowerStartMoliereRadius(-10.f),
    m_middlePostShowerStartShowerStartMoliereRadius(-10.f),
    m_maxPostShowerStartShowerStartMoliereRadius(-10.f),
    m_positiveOpeningAngleW(-10.f),
    m_negativeOpeningAngleW(-10.f),
    m_maxOpeningAngleW(-10.f),
    m_showerApexLW(-9999.f),
    m_minShowerApexL(-9999.f),
    m_middleShowerApexL(-9999.f),
    m_maxShowerApexL(-9999.f),
    m_showerApexTW(-10.f),
    m_minShowerApexT(-10.f),
    m_middleShowerApexT(-10.f),
    m_maxShowerApexT(-10.f),
    m_foundHitRatioW(-0.5f),
    m_minFoundHitRatio(-0.5f),
    m_middleFoundHitRatio(-0.5f),
    m_maxFoundHitRatio(-0.5f),
    m_fitShowerStartLW(-10.f),
    m_fitShowerStartTW(-10.f),
    m_postShowerStartMinHalfOpeningAngle(-10.f),
    m_postShowerStartMaxHalfOpeningAngle(-10.f),
    m_postShowerStartOpeningAngle(-10.f),
    m_postShowerStartOpeningAngleAsymmetry(-10.f),
    m_postShowerStartMeanTransverseAngle(-10.f),
    m_postShowerStartMeanLWeightedTransverseAngle(-10.f),
    m_postShowerStartMeanRadialDistance(-10.f),
    m_postShowerStartRadialDistanceSigma(-10.f),
    m_postShowerStartEnergyWeightedMeanRadialDistance(-10.f),
    m_postShowerStartEstimatedMoliereRadius(-10.f),
    m_postShowerStartLWeightedMeanRadialDistance(-10.f),
    m_postShowerStartLWeightedRadialDistanceSigma(-10.f),
    m_postShowerStartInitialGapSize(-10.f),
    m_postShowerStartMaxGapSize(-10.f),
    m_initialRegionDistanceToNuVertex(-10.f),
    m_initialRegionDistanceInGaps(-10.f),
    m_initialRegionMaxGapSize(-10.f),
    m_initialGapSizeW(-10.f),
    m_minInitialGapSize(-10.f),
    m_middleInitialGapSize(-10.f),
    m_maxInitialGapSize(-10.f),
    m_largestGapSizeW(-10.f),
    m_minLargestGapSize(-10.f),
    m_middleLargestGapSize(-10.f),
    m_maxLargestGapSize(-10.f),
    m_largestProjectedGapSizeW(-10.f),
    m_minLargestProjectedGapSize(-10.f),
    m_middleLargestProjectedGapSize(-10.f),
    m_maxLargestProjectedGapSize(-10.f),
    m_hitLineDensityW(-10.f),
    m_minHitLineDensity(-10.f),
    m_middleHitLineDensity(-10.f),
    m_maxHitLineDensity(-10.f),
    m_nViewsWithAmbiguousHits(-10),
    m_nAmbiguousHits2D(-10),
    m_minNAmbiguousHits(-10),
    m_maxNAmbiguousHits(-10),
    m_ambiguousHitUnaccountedEnergyU(-10.f),
    m_ambiguousHitUnaccountedEnergyV(-10.f),
    m_ambiguousHitUnaccountedEnergyW(-10.f),
    m_ambiguousHitMinUnaccountedEnergy(-10.f),
    m_ambiguousHitMaxUnaccountedEnergy(-10.f),
    m_ambiguousHitShowerEnergyRatioU(-10.f),
    m_ambiguousHitShowerEnergyRatioV(-10.f),
    m_ambiguousHitShowerEnergyRatioW(-10.f),
    m_ambiguousHitMinShowerEnergyRatio(-10.f),
    m_ambiguousHitMaxShowerEnergyRatio(-10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArConnectionPathwayHelper::ElectronTreeVariables &LArConnectionPathwayHelper::ElectronTreeVariables::operator=(const ElectronTreeVariables &rhs)
{
    m_nConnectionPathways = rhs.m_nConnectionPathways;
    m_showerStartX = rhs.m_showerStartX;
    m_showerStartY = rhs.m_showerStartY;
    m_showerStartZ = rhs.m_showerStartZ;
    m_pathwayLengthMin = rhs.m_pathwayLengthMin;
    m_pathwayLengthMiddle = rhs.m_pathwayLengthMiddle;
    m_pathwayLengthMax = rhs.m_pathwayLengthMax;
    m_pathwayShowerStartDelta = rhs.m_pathwayShowerStartDelta;
    m_pathwayMaxScatteringAngle = rhs.m_pathwayMaxScatteringAngle;
    m_minShowerStartPathwayScatteringAngle2D = rhs.m_minShowerStartPathwayScatteringAngle2D;
    m_middleShowerStartPathwayScatteringAngle2D = rhs.m_middleShowerStartPathwayScatteringAngle2D;
    m_maxShowerStartPathwayScatteringAngle2D = rhs.m_maxShowerStartPathwayScatteringAngle2D;
    m_pathwayEnergyMeanU = rhs.m_pathwayEnergyMeanU;
    m_pathwayEnergyMeanV = rhs.m_pathwayEnergyMeanV;
    m_pathwayEnergyMeanW = rhs.m_pathwayEnergyMeanW;
    m_pathwayEnergySigmaU = rhs.m_pathwayEnergySigmaU;
    m_pathwayEnergySigmaV = rhs.m_pathwayEnergySigmaV;
    m_pathwayEnergySigmaW = rhs.m_pathwayEnergySigmaW;
    m_pathwayMinEnergyMean = rhs.m_pathwayMinEnergyMean;
    m_pathwayMiddleEnergyMean = rhs.m_pathwayMiddleEnergyMean;
    m_pathwayMaxEnergyMean = rhs.m_pathwayMaxEnergyMean;
    m_pathwayMinEnergySigma = rhs.m_pathwayMinEnergySigma;
    m_pathwayMiddleEnergySigma = rhs.m_pathwayMiddleEnergySigma;
    m_pathwayMaxEnergySigma = rhs.m_pathwayMaxEnergySigma;
    m_postShowerStartLength = rhs.m_postShowerStartLength;
    m_postShowerStartNHits = rhs.m_postShowerStartNHits;
    m_postShowerStartNHitsW = rhs.m_postShowerStartNHitsW;
    m_minNPostShowerStartHits = rhs.m_minNPostShowerStartHits;
    m_middleNPostShowerStartHits = rhs.m_middleNPostShowerStartHits;
    m_maxNPostShowerStartHits = rhs.m_maxNPostShowerStartHits;
    m_postShowerStartScatterAngle = rhs.m_postShowerStartScatterAngle;
    m_postShowerStartScatterAngleW = rhs.m_postShowerStartScatterAngleW;
    m_minPostShowerStartScatterAngle = rhs.m_minPostShowerStartScatterAngle;
    m_middlePostShowerStartScatterAngle = rhs.m_middlePostShowerStartScatterAngle;
    m_maxPostShowerStartScatterAngle = rhs.m_maxPostShowerStartScatterAngle;
    m_postShowerStartOpeningAngleW = rhs.m_postShowerStartOpeningAngleW;
    m_minPostShowerStartOpeningAngle = rhs.m_minPostShowerStartOpeningAngle;
    m_middlePostShowerStartOpeningAngle = rhs.m_middlePostShowerStartOpeningAngle;
    m_maxPostShowerStartOpeningAngle = rhs.m_maxPostShowerStartOpeningAngle;
    m_postShowerStartOpeningAngleAsymmetryW = rhs.m_postShowerStartOpeningAngleAsymmetryW;
    m_minPostShowerStartOpeningAngleAsymmetry = rhs.m_minPostShowerStartOpeningAngleAsymmetry;
    m_middlePostShowerStartOpeningAngleAsymmetry = rhs.m_middlePostShowerStartOpeningAngleAsymmetry;
    m_maxPostShowerStartOpeningAngleAsymmetry = rhs.m_maxPostShowerStartOpeningAngleAsymmetry;
    m_postShowerStartNuVertexHitAsymmetryW = rhs.m_postShowerStartNuVertexHitAsymmetryW;
    m_minPostShowerStartNuVertexHitAsymmetry = rhs.m_minPostShowerStartNuVertexHitAsymmetry;
    m_middlePostShowerStartNuVertexHitAsymmetry = rhs.m_middlePostShowerStartNuVertexHitAsymmetry;
    m_maxPostShowerStartNuVertexHitAsymmetry = rhs.m_maxPostShowerStartNuVertexHitAsymmetry;
    m_postShowerStartNuVertexEnergyAsymmetryW = rhs.m_postShowerStartNuVertexEnergyAsymmetryW;
    m_minPostShowerStartNuVertexEnergyAsymmetry = rhs.m_minPostShowerStartNuVertexEnergyAsymmetry;
    m_middlePostShowerStartNuVertexEnergyAsymmetry = rhs.m_middlePostShowerStartNuVertexEnergyAsymmetry;
    m_maxPostShowerStartNuVertexEnergyAsymmetry = rhs.m_maxPostShowerStartNuVertexEnergyAsymmetry;
    m_postShowerStartShowerStartHitAsymmetryW = rhs.m_postShowerStartShowerStartHitAsymmetryW;
    m_minPostShowerStartShowerStartHitAsymmetry = rhs.m_minPostShowerStartShowerStartHitAsymmetry;
    m_middlePostShowerStartShowerStartHitAsymmetry = rhs.m_middlePostShowerStartShowerStartHitAsymmetry;
    m_maxPostShowerStartShowerStartHitAsymmetry = rhs.m_maxPostShowerStartShowerStartHitAsymmetry;
    m_postShowerStartShowerStartEnergyAsymmetryW = rhs.m_postShowerStartShowerStartEnergyAsymmetryW;
    m_minPostShowerStartShowerStartEnergyAsymmetry = rhs.m_minPostShowerStartShowerStartEnergyAsymmetry;
    m_middlePostShowerStartShowerStartEnergyAsymmetry = rhs.m_middlePostShowerStartShowerStartEnergyAsymmetry;
    m_maxPostShowerStartShowerStartEnergyAsymmetry = rhs.m_maxPostShowerStartShowerStartEnergyAsymmetry;        
    m_postShowerStartNuVertexMeanRadialDistanceW = rhs.m_postShowerStartNuVertexMeanRadialDistanceW;
    m_minPostShowerStartNuVertexMeanRadialDistance = rhs.m_minPostShowerStartNuVertexMeanRadialDistance;
    m_middlePostShowerStartNuVertexMeanRadialDistance = rhs.m_middlePostShowerStartNuVertexMeanRadialDistance;
    m_maxPostShowerStartNuVertexMeanRadialDistance = rhs.m_maxPostShowerStartNuVertexMeanRadialDistance;
    m_postShowerStartNuVertexEnergyWeightedMeanRadialDistanceW = rhs.m_postShowerStartNuVertexEnergyWeightedMeanRadialDistanceW;
    m_minPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = rhs.m_minPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    m_middlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance = rhs.m_middlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = rhs.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    m_postShowerStartShowerStartMeanRadialDistanceW = rhs.m_postShowerStartShowerStartMeanRadialDistanceW;
    m_minPostShowerStartShowerStartMeanRadialDistance = rhs.m_minPostShowerStartShowerStartMeanRadialDistance;
    m_middlePostShowerStartShowerStartMeanRadialDistance = rhs.m_middlePostShowerStartShowerStartMeanRadialDistance;
    m_maxPostShowerStartShowerStartMeanRadialDistance = rhs.m_maxPostShowerStartShowerStartMeanRadialDistance;
    m_postShowerStartShowerStartEnergyWeightedMeanRadialDistanceW = rhs.m_postShowerStartShowerStartEnergyWeightedMeanRadialDistanceW;
    m_minPostShowerStartShowerStartEnergyWeightedMeanRadialDistance = rhs.m_minPostShowerStartShowerStartEnergyWeightedMeanRadialDistance;
    m_middlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance = rhs.m_middlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance;
    m_maxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance = rhs.m_maxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance;
    m_postShowerStartNuVertexMoliereRadiusW = rhs.m_postShowerStartNuVertexMoliereRadiusW;
    m_minPostShowerStartNuVertexMoliereRadius = rhs.m_minPostShowerStartNuVertexMoliereRadius;
    m_middlePostShowerStartNuVertexMoliereRadius = rhs.m_middlePostShowerStartNuVertexMoliereRadius;
    m_maxPostShowerStartNuVertexMoliereRadius = rhs.m_maxPostShowerStartNuVertexMoliereRadius;
    m_postShowerStartShowerStartMoliereRadiusW = rhs.m_postShowerStartShowerStartMoliereRadiusW;
    m_minPostShowerStartShowerStartMoliereRadius = rhs.m_minPostShowerStartShowerStartMoliereRadius;
    m_middlePostShowerStartShowerStartMoliereRadius = rhs.m_middlePostShowerStartShowerStartMoliereRadius;
    m_maxPostShowerStartShowerStartMoliereRadius = rhs.m_maxPostShowerStartShowerStartMoliereRadius;
    m_positiveOpeningAngleW = rhs.m_positiveOpeningAngleW;
    m_negativeOpeningAngleW = rhs.m_negativeOpeningAngleW;
    m_maxOpeningAngleW = rhs.m_maxOpeningAngleW;
    m_showerApexLW = rhs.m_showerApexLW;
    m_minShowerApexL = rhs.m_minShowerApexL;
    m_middleShowerApexL = rhs.m_middleShowerApexL;
    m_maxShowerApexL = rhs.m_maxShowerApexL;
    m_showerApexTW = rhs.m_showerApexTW;
    m_minShowerApexT = rhs.m_minShowerApexT;
    m_middleShowerApexT = rhs.m_middleShowerApexT;
    m_maxShowerApexT = rhs.m_maxShowerApexT;
    m_foundHitRatioW = rhs.m_foundHitRatioW;
    m_minFoundHitRatio = rhs.m_minFoundHitRatio;
    m_middleFoundHitRatio = rhs.m_middleFoundHitRatio;
    m_maxFoundHitRatio = rhs.m_maxFoundHitRatio;
    m_fitShowerStartLW = rhs.m_fitShowerStartLW;
    m_fitShowerStartTW = rhs.m_fitShowerStartTW;
    m_postShowerStartMinHalfOpeningAngle = rhs.m_postShowerStartMinHalfOpeningAngle;
    m_postShowerStartMaxHalfOpeningAngle = rhs.m_postShowerStartMaxHalfOpeningAngle;
    m_postShowerStartOpeningAngle = rhs.m_postShowerStartOpeningAngle;
    m_postShowerStartOpeningAngleAsymmetry = rhs.m_postShowerStartOpeningAngleAsymmetry;
    m_postShowerStartMeanTransverseAngle = rhs.m_postShowerStartMeanTransverseAngle;
    m_postShowerStartMeanLWeightedTransverseAngle = rhs.m_postShowerStartMeanLWeightedTransverseAngle;
    m_postShowerStartMeanRadialDistance = rhs.m_postShowerStartMeanRadialDistance;
    m_postShowerStartRadialDistanceSigma = rhs.m_postShowerStartRadialDistanceSigma;
    m_postShowerStartEnergyWeightedMeanRadialDistance = rhs.m_postShowerStartEnergyWeightedMeanRadialDistance;
    m_postShowerStartEstimatedMoliereRadius = rhs.m_postShowerStartEstimatedMoliereRadius;
    m_postShowerStartLWeightedMeanRadialDistance = rhs.m_postShowerStartLWeightedMeanRadialDistance;
    m_postShowerStartLWeightedRadialDistanceSigma = rhs.m_postShowerStartLWeightedRadialDistanceSigma;
    m_postShowerStartInitialGapSize = rhs.m_postShowerStartInitialGapSize;
    m_postShowerStartMaxGapSize = rhs.m_postShowerStartMaxGapSize;
    m_initialRegionDistanceToNuVertex = rhs.m_initialRegionDistanceToNuVertex;
    m_initialRegionDistanceInGaps = rhs.m_initialRegionDistanceInGaps;
    m_initialRegionMaxGapSize = rhs.m_initialRegionMaxGapSize;
    m_initialGapSizeW = rhs.m_initialGapSizeW;
    m_minInitialGapSize = rhs.m_minInitialGapSize;
    m_middleInitialGapSize = rhs.m_middleInitialGapSize;
    m_maxInitialGapSize = rhs.m_maxInitialGapSize;
    m_largestGapSizeW = rhs.m_largestGapSizeW;
    m_minLargestGapSize = rhs.m_minLargestGapSize;
    m_middleLargestGapSize = rhs.m_middleLargestGapSize;
    m_maxLargestGapSize = rhs.m_maxLargestGapSize;
    m_largestProjectedGapSizeW = rhs.m_largestProjectedGapSizeW;
    m_minLargestProjectedGapSize = rhs.m_minLargestProjectedGapSize;
    m_middleLargestProjectedGapSize = rhs.m_middleLargestProjectedGapSize;
    m_maxLargestProjectedGapSize = rhs.m_maxLargestProjectedGapSize;
    m_hitLineDensityW = rhs.m_hitLineDensityW;
    m_minHitLineDensity = rhs.m_minHitLineDensity;
    m_middleHitLineDensity = rhs.m_middleHitLineDensity;
    m_maxHitLineDensity = rhs.m_maxHitLineDensity;
    m_nViewsWithAmbiguousHits = rhs.m_nViewsWithAmbiguousHits;
    m_nAmbiguousHits2D = rhs.m_nAmbiguousHits2D;
    m_minNAmbiguousHits = rhs.m_minNAmbiguousHits;
    m_maxNAmbiguousHits = rhs.m_maxNAmbiguousHits;
    m_ambiguousHitUnaccountedEnergyU = rhs.m_ambiguousHitUnaccountedEnergyU;
    m_ambiguousHitUnaccountedEnergyV = rhs.m_ambiguousHitUnaccountedEnergyV;
    m_ambiguousHitUnaccountedEnergyW = rhs.m_ambiguousHitUnaccountedEnergyW;
    m_ambiguousHitMinUnaccountedEnergy = rhs.m_ambiguousHitMinUnaccountedEnergy;
    m_ambiguousHitMaxUnaccountedEnergy = rhs.m_ambiguousHitMaxUnaccountedEnergy;
    m_ambiguousHitShowerEnergyRatioU = rhs.m_ambiguousHitShowerEnergyRatioU;
    m_ambiguousHitShowerEnergyRatioV = rhs.m_ambiguousHitShowerEnergyRatioV;
    m_ambiguousHitShowerEnergyRatioW = rhs.m_ambiguousHitShowerEnergyRatioW;
    m_ambiguousHitMinShowerEnergyRatio = rhs.m_ambiguousHitMinShowerEnergyRatio;
    m_ambiguousHitMaxShowerEnergyRatio = rhs.m_ambiguousHitMaxShowerEnergyRatio;

    return *this;
}       

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::FillElectronTreeVariables(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW, const CartesianVector &nuVertexPosition, 
    const CaloHitList *const pCaloHitListU, const CaloHitList *const pCaloHitListV, const CaloHitList *const pCaloHitListW,
    const Consistency &consistency, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    CartesianVector minShowerStart3D(0.f, 0.f, 0.f), middleShowerStart3D(0.f, 0.f, 0.f), maxShowerStart3D(0.f, 0.f, 0.f);

    ProtoShowerMatch protoShowerMatch(protoShowerU, protoShowerV, protoShowerW, consistency);
    CartesianPointVector showerStarts3D;

    if (!LArConnectionPathwayHelper::FindShowerStarts3D(pAlgorithm, pShowerPfo, protoShowerMatch, nuVertexPosition, 
        showerStarts3D))
    {
        /////////////////////////////////
        std::cout << "COULD NOT FIND ANY SHOWER VERTICES!" << std::endl;
        /////////////////////////////////
        return;
    }

    minShowerStart3D = showerStarts3D.at(0);
    middleShowerStart3D = showerStarts3D.at(1);
    maxShowerStart3D = showerStarts3D.at(2);

    electronTreeVariables.m_showerStartX = minShowerStart3D.GetX();
    electronTreeVariables.m_showerStartY = minShowerStart3D.GetY();
    electronTreeVariables.m_showerStartZ = minShowerStart3D.GetZ();

    std::cout << "electronTreeVariables.m_showerStartX: " << electronTreeVariables.m_showerStartX << std::endl;
    std::cout << "electronTreeVariables.m_showerStartY: " << electronTreeVariables.m_showerStartY << std::endl;
    std::cout << "electronTreeVariables.m_showerStartZ: " << electronTreeVariables.m_showerStartZ << std::endl;

    /////////////////////////////////
    /*
    std::cout << "(minShowerStart3D - nuVertexPosition).GetMagnitude(): " << (minShowerStart3D - nuVertexPosition).GetMagnitude() << std::endl;
    std::cout << "(middleShowerStart3D - nuVertexPosition).GetMagnitude(): " << (middleShowerStart3D - nuVertexPosition).GetMagnitude() << std::endl;
    std::cout << "(maxShowerStart3D - nuVertexPosition).GetMagnitude(): " << (maxShowerStart3D - nuVertexPosition).GetMagnitude() << std::endl;
    
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &minShowerStart3D, "MIN", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &middleShowerStart3D, "MIDDLE", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &maxShowerStart3D, "MAX", BLUE, 2);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "U SPINE", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "V SPINE", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "W SPINE", GREEN);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    //std::cout << "//////////////////////////////" << std::endl;
    //std::cout << "NOW FILLING CONNECTION PATHWAY VARIABLES..." << std::endl;
    LArConnectionPathwayHelper::FillConnectionPathwayVariables(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, pShowerPfo, nuVertexPosition, minShowerStart3D,
        middleShowerStart3D, maxShowerStart3D, electronTreeVariables);
    //min, middle, max


    //std::cout << "//////////////////////////////" << std::endl;
    //std::cout << "NOW FILLING POST SHOWER START VARIABLES..." << std::endl;
    LArConnectionPathwayHelper::FillPostShowerStartVariables(pAlgorithm, pShowerPfo, nuVertexPosition, protoShowerU, protoShowerV, protoShowerW, minShowerStart3D, electronTreeVariables);
    // middle


    //std::cout << "//////////////////////////////" << std::endl;
    //std::cout << "NOW FILLING INITIAL REGION VARIABLES..." << std::endl;
    LArConnectionPathwayHelper::FillInitialRegionVariables(pAlgorithm, nuVertexPosition, protoShowerU, protoShowerV, protoShowerW, 
        minShowerStart3D, electronTreeVariables);
    // middle


    //std::cout << "//////////////////////////////" << std::endl;
    //std::cout << "NOW FILLING AMBIGUOUS HIT VARIABLES..." << std::endl;
    LArConnectionPathwayHelper::FillAmbiguousHitVariables(pAlgorithm, pShowerPfo, protoShowerU, protoShowerV, protoShowerW, nuVertexPosition, minShowerStart3D, pCaloHitListU, 
        pCaloHitListV, pCaloHitListW, electronTreeVariables);
    //middle
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const CartesianVector &nuVertexPosition,
    CartesianVector &middleShowerStart3D, const CaloHitList *const pCaloHitListU, const CaloHitList *const pCaloHitListV, const CaloHitList *const pCaloHitListW,
    LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    int nViewsWithAmbiguousHits(0);

    const int nAmbiguousHitsU(protoShowerU.m_ambiguousHitList.size());
    nViewsWithAmbiguousHits += (nAmbiguousHitsU == 0) ? 0 : 1;

    const int nAmbiguousHitsV(protoShowerV.m_ambiguousHitList.size());
    nViewsWithAmbiguousHits += (nAmbiguousHitsV == 0) ? 0 : 1;

    const int nAmbiguousHitsW(protoShowerW.m_ambiguousHitList.size());
    nViewsWithAmbiguousHits += (nAmbiguousHitsW == 0) ? 0 : 1;

    electronTreeVariables.m_nViewsWithAmbiguousHits = nViewsWithAmbiguousHits;
    //std::cout << "electronTreeVariables.m_nViewsWithAmbiguousHits: " << electronTreeVariables.m_nViewsWithAmbiguousHits << std::endl;

    const int nAmbiguousHits2D(nAmbiguousHitsU + nAmbiguousHitsV + nAmbiguousHitsW);
    electronTreeVariables.m_nAmbiguousHits2D = nAmbiguousHits2D;
    //std::cout << "electronTreeVariables.m_nAmbiguousHits2D: " << electronTreeVariables.m_nAmbiguousHits2D << std::endl;

    const int minNAmbiguousHits(std::min(std::min(nAmbiguousHitsU, nAmbiguousHitsV), nAmbiguousHitsW));
    electronTreeVariables.m_minNAmbiguousHits = minNAmbiguousHits;
    //std::cout << "electronTreeVariables.m_minNAmbiguousHits: " << electronTreeVariables.m_minNAmbiguousHits << std::endl;

    const int maxNAmbiguousHits(std::max(std::max(nAmbiguousHitsU, nAmbiguousHitsV), nAmbiguousHitsW));
    electronTreeVariables.m_maxNAmbiguousHits = maxNAmbiguousHits;
    //std::cout << "electronTreeVariables.m_maxNAmbiguousHits: " << electronTreeVariables.m_maxNAmbiguousHits << std::endl;

    //////////////////////////////
    /*
    for (const CaloHit *const pCaloHit : protoShowerU.m_ambiguousHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ambiguous hit U", GRAY, 2);
    }

    for (const CaloHit *const pCaloHit : protoShowerV.m_ambiguousHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ambiguous hit V", GRAY, 2);
    }

    for (const CaloHit *const pCaloHit : protoShowerW.m_ambiguousHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ambiguous hit W", GRAY, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////////

    if ((nAmbiguousHitsU == 0) && (nAmbiguousHitsV == 0) && (nAmbiguousHitsW == 0))
    {
        //std::cout << "no ambiguous hits, returning..." << std::endl;
        return false;
    }

    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    try
    {
        // try 10 so i can make it consistent?
        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const CartesianVector &startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector &startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector &startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
            protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
        {
            if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
                protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, protoShowerU.m_connectionPathway.m_startDirection, protoShowerV.m_connectionPathway.m_startDirection, 
                protoShowerW.m_connectionPathway.m_startDirection, connectionPathwayDirection3D))
            {
                //std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
                return false;
            }
        }

        bool found(false);
        float minUnaccountedEnergy(std::numeric_limits<float>::max()), maxUnaccountedEnergy(-std::numeric_limits<float>::max());
        float minShowerEnergyRatio(std::numeric_limits<float>::max()), maxShowerEnergyRatio(-std::numeric_limits<float>::max());

        float unaccountedHitEnergyU(-10.f), showerEnergyRatioU(-10.f);
        if (LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pAlgorithm, protoShowerU, TPC_VIEW_U, nuVertexPosition, connectionPathwayDirection3D, 
            pCaloHitListU, unaccountedHitEnergyU, showerEnergyRatioU))
        {
            found = true;
            minUnaccountedEnergy = std::min(minUnaccountedEnergy, unaccountedHitEnergyU);
            maxUnaccountedEnergy = std::max(maxUnaccountedEnergy, unaccountedHitEnergyU);
            minShowerEnergyRatio = std::min(minShowerEnergyRatio, showerEnergyRatioU);
            maxShowerEnergyRatio = std::max(maxShowerEnergyRatio, showerEnergyRatioU);
        }

        float unaccountedHitEnergyV(-10.f), showerEnergyRatioV(-10.f);
        if (LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pAlgorithm, protoShowerV, TPC_VIEW_V, nuVertexPosition, connectionPathwayDirection3D, 
            pCaloHitListV, unaccountedHitEnergyV, showerEnergyRatioV))
        {
            found = true;
            minUnaccountedEnergy = std::min(minUnaccountedEnergy, unaccountedHitEnergyV);
            maxUnaccountedEnergy = std::max(maxUnaccountedEnergy, unaccountedHitEnergyV);
            minShowerEnergyRatio = std::min(minShowerEnergyRatio, showerEnergyRatioV);
            maxShowerEnergyRatio = std::max(maxShowerEnergyRatio, showerEnergyRatioV);
        }

        float unaccountedHitEnergyW(-10.f), showerEnergyRatioW(-10.f);
        if (LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pAlgorithm, protoShowerW, TPC_VIEW_W, nuVertexPosition, connectionPathwayDirection3D, 
            pCaloHitListW, unaccountedHitEnergyW, showerEnergyRatioW))
        {
            found = true;
            minUnaccountedEnergy = std::min(minUnaccountedEnergy, unaccountedHitEnergyW);
            maxUnaccountedEnergy = std::max(maxUnaccountedEnergy, unaccountedHitEnergyW);
            minShowerEnergyRatio = std::min(minShowerEnergyRatio, showerEnergyRatioW);
            maxShowerEnergyRatio = std::max(maxShowerEnergyRatio, showerEnergyRatioW);
        }

        if (!found)
        {
            minUnaccountedEnergy = -10.f;
            maxUnaccountedEnergy = -10.f;
            minShowerEnergyRatio = -10.f;
            maxShowerEnergyRatio = -10.f;
        }

        electronTreeVariables.m_ambiguousHitUnaccountedEnergyU = unaccountedHitEnergyU;
        //std::cout << "electronTreeVariables.m_ambiguousHitUnaccountedEnergyU: " << electronTreeVariables.m_ambiguousHitUnaccountedEnergyU << std::endl;

        electronTreeVariables.m_ambiguousHitUnaccountedEnergyV = unaccountedHitEnergyV;
        //std::cout << "electronTreeVariables.m_ambiguousHitUnaccountedEnergyV: " << electronTreeVariables.m_ambiguousHitUnaccountedEnergyV << std::endl;

        electronTreeVariables.m_ambiguousHitUnaccountedEnergyW = unaccountedHitEnergyW;
        //std::cout << "electronTreeVariables.m_ambiguousHitUnaccountedEnergyW: " << electronTreeVariables.m_ambiguousHitUnaccountedEnergyW << std::endl;

        electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy = minUnaccountedEnergy;
        //std::cout << "electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy: " << electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy << std::endl;

        electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy = maxUnaccountedEnergy;
        //std::cout << "electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy: " << electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy << std::endl;

        electronTreeVariables.m_ambiguousHitShowerEnergyRatioU = showerEnergyRatioU;
        //std::cout << "electronTreeVariables.m_ambiguousHitShowerEnergyRatioU: " << electronTreeVariables.m_ambiguousHitShowerEnergyRatioU << std::endl;

        electronTreeVariables.m_ambiguousHitShowerEnergyRatioV = showerEnergyRatioV;
        //std::cout << "electronTreeVariables.m_ambiguousHitShowerEnergyRatioV: " << electronTreeVariables.m_ambiguousHitShowerEnergyRatioV << std::endl;

        electronTreeVariables.m_ambiguousHitShowerEnergyRatioW = showerEnergyRatioW;
        //std::cout << "electronTreeVariables.m_ambiguousHitShowerEnergyRatioW: " << electronTreeVariables.m_ambiguousHitShowerEnergyRatioW << std::endl;

        electronTreeVariables.m_ambiguousHitMinShowerEnergyRatio = minShowerEnergyRatio;
        //std::cout << "electronTreeVariables.m_ambiguousHitMinShowerEnergyRatio: " << electronTreeVariables.m_ambiguousHitMinShowerEnergyRatio << std::endl;

        electronTreeVariables.m_ambiguousHitMaxShowerEnergyRatio = maxShowerEnergyRatio;
        //std::cout << "electronTreeVariables.m_ambiguousHitMaxShowerEnergyRatio: " << electronTreeVariables.m_ambiguousHitMaxShowerEnergyRatio << std::endl;
    }
    catch (...)
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const HitType hitType, const CartesianVector &nuVertexPosition, const CartesianVector &connectionPathwayDirection3D,
    const CaloHitList *const pCaloHitList, float &unaccountedHitEnergy, float &showerEnergyRatio)
{
    const CartesianVector projectedNuVertex(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, hitType));

    // Find the owner connection pathways of the ambiguous hits
    std::map<int, CaloHitList> ambiguousHitSpinesTemp;
    CaloHitList hitsToExcludeInEnergyCalcs; // to avoid double counting the energy

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (std::find(protoShower.m_ambiguousHitList.begin(), protoShower.m_ambiguousHitList.end(), pCaloHit) != protoShower.m_ambiguousHitList.end())
            continue;

        int count(0);

        // A hit can be in more than one spine
        for (unsigned int i = 0; i < protoShower.m_ambiguousDirectionVector.size(); ++i)
        {
            const CartesianVector &significantDirection(protoShower.m_ambiguousDirectionVector[i]);
            const CartesianVector displacement(pCaloHit->GetPositionVector() - projectedNuVertex);
            const float thisT(significantDirection.GetCrossProduct(displacement).GetMagnitude());
            const float thisL(significantDirection.GetDotProduct(displacement));

            if ((thisL > 0.f) && (thisT < 0.75f))
            {
                ++count;
                ambiguousHitSpinesTemp[i].push_back(pCaloHit);
            }

            if (count == 2)
                hitsToExcludeInEnergyCalcs.push_back(pCaloHit);
        }
    }

    if (ambiguousHitSpinesTemp.empty())
    {
        //std::cout << "NO AMBIGUOUS SPINES FOUND" << std::endl;
        return false;
    }

    /////////////////////////////////
    /*
    for (const auto &entry : ambiguousHitSpinesTemp)
    {
        for (const CaloHit *const pCaloHit : entry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ambiguous spine hit", DARKGREEN, 2);
        }
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    ///////////////////////////////// 
 
    // Find continuous pathways i.e. make sure other connection pathway objects are decent
    std::map<int, CaloHitList> ambiguousHitSpines;

    for (const auto &entry : ambiguousHitSpinesTemp)
    {
        CaloHitList continuousSpine(LArConnectionPathwayHelper::FindAmbiguousContinuousSpine(entry.second, protoShower.m_ambiguousHitList, projectedNuVertex));

        if (continuousSpine.size() > 0)
            ambiguousHitSpines[entry.first] = continuousSpine;
    }

    if (ambiguousHitSpines.empty())
    {
        //std::cout << "couldn't find a connected pathway" << std::endl;
        return false;
    }

    /////////////////////////////////
    /*
    for (const auto &entry : ambiguousHitSpines)
    {
        for (const CaloHit *const pCaloHit : entry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "connected ambiguous spine hit", TEAL, 2);
        }
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    // We can't do a de/dx without a 3D direction :( 
    // Which is really hard to obtain for the other particles in the event... so we might have to do mean energy? 
    // But we can make sure that we only consider a certain length of true spine (lets say 5cm)
    // Pick this to avoid going into the shower i.e. hopefully shower dedx will be a bit constant

    const CartesianVector seedPoint3D(nuVertexPosition + (connectionPathwayDirection3D * 3.f));
    const CartesianVector seedPoint(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), seedPoint3D, hitType));
    const CartesianVector startDirection((seedPoint - projectedNuVertex).GetUnitVector());
    const float lRange((seedPoint - projectedNuVertex).GetMagnitude());

    // Need to find max start position...
    float ambiguousHitEnergyMean(0.f);
    float startL(-std::numeric_limits<float>::max());

    for (const CaloHit *const pAmbiguousCaloHit : protoShower.m_ambiguousHitList)
    {
        const float thisT(startDirection.GetCrossProduct(pAmbiguousCaloHit->GetPositionVector() - projectedNuVertex).GetMagnitude());
        const float thisL(startDirection.GetDotProduct(pAmbiguousCaloHit->GetPositionVector() - projectedNuVertex));

        if ((thisL > startL) && (thisT < 0.75f))
            startL = thisL;

        ambiguousHitEnergyMean += pAmbiguousCaloHit->GetElectromagneticEnergy() * 1000.f;
    }

    if (startL < 0.f)
    {
        //std::cout << "(startL < 0.f) so leaving..." << std::endl;
        return false;
    }

    ambiguousHitEnergyMean /= protoShower.m_ambiguousHitList.size();

    // Get mean energy of other pathways, avoiding the double counting hits (ambiguous to the ambiguous hits)
    float otherEnergyMeanSum(0.f);

    for (const auto &entry : ambiguousHitSpines)
    {
        int nOtherEnergyHits(0);
        float otherEnergyMean(0.f);

        for (const CaloHit *const pOtherCaloHit : entry.second)
        {
            if (std::find(hitsToExcludeInEnergyCalcs.begin(), hitsToExcludeInEnergyCalcs.end(), pOtherCaloHit) != hitsToExcludeInEnergyCalcs.end())
                continue;

            otherEnergyMean += pOtherCaloHit->GetElectromagneticEnergy() * 1000.f;
            ++nOtherEnergyHits;
        }

        if (nOtherEnergyHits == 0)
            continue;

        otherEnergyMean /= static_cast<float>(nOtherEnergyHits);
        otherEnergyMeanSum += otherEnergyMean;
    }

    // Get the spine mean energy

    float spineEnergyMean(0.f);
    int nSpineEnergyHits(0);

    for (const CaloHit *const pSpineCaloHit : protoShower.m_spineHitList)
    {
        if (std::find(protoShower.m_ambiguousHitList.begin(), protoShower.m_ambiguousHitList.end(), pSpineCaloHit) != protoShower.m_ambiguousHitList.end())
            continue;

        const float thisL(startDirection.GetDotProduct(pSpineCaloHit->GetPositionVector() - projectedNuVertex));

        if ((thisL > startL) && (thisL < (startL + lRange)))
        {
            spineEnergyMean += pSpineCaloHit->GetElectromagneticEnergy() * 1000.f;
            ++nSpineEnergyHits;
        }

        ///////////////////////////
        /*
        const CartesianVector &hitPosition(pSpineCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "energy calc spine hit", DARKYELLOW, 2);
        */
        ///////////////////////////
    }

    if (nSpineEnergyHits == 0)
    {
        //std::cout << "no spine hits to use in energy calculation" << std::endl;
        return false;
    }

    ///////////////////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    ///////////////////////////////

    spineEnergyMean /= static_cast<float>(nSpineEnergyHits);

    unaccountedHitEnergy = ambiguousHitEnergyMean - otherEnergyMeanSum - spineEnergyMean;
    showerEnergyRatio = (ambiguousHitEnergyMean - otherEnergyMeanSum) / spineEnergyMean;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArConnectionPathwayHelper::FindAmbiguousContinuousSpine(const CaloHitList &caloHitList, const CaloHitList &ambiguousHitList, 
    const CartesianVector &projectedNuVertexPosition)
{
    CaloHitList continuousHitList;

    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(projectedNuVertexPosition));

    for (unsigned int i = 0; i < caloHitVector.size(); ++i)
    {
        CaloHitList connectedHitList;
        connectedHitList.push_back(caloHitVector[i]);

        if (LArClusterHelper::GetClosestDistance(connectedHitList.front()->GetPositionVector(), ambiguousHitList) > 1.f)
            continue;

        bool found(true);

        while(found)
        {
            found = false;

            for (unsigned int j = (i + 1); j < caloHitVector.size(); ++j)
            {
                const CaloHit *const pCaloHit(caloHitVector[j]);

                if (std::find(connectedHitList.begin(), connectedHitList.end(), pCaloHit) != connectedHitList.end())
                    continue;

                if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), connectedHitList) < 1.f)
                {
                    // to avoid ends of tracks
                    if (static_cast<float>(connectedHitList.size()) < static_cast<float>(caloHitVector.size() * 0.8f))
                    {
                        found = true;
                        connectedHitList.push_back(pCaloHit);
                    }

                    break;
                }
            }
        }

        if (connectedHitList.size() >= 2)
        {
            continuousHitList.insert(continuousHitList.begin(), connectedHitList.begin(), connectedHitList.end());
            break;
        }
    }

    return continuousHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ALSO ADD DISTANCE FROM NU VERTEX

bool LArConnectionPathwayHelper::FillInitialRegionVariables(pandora::Algorithm *const pAlgorithm, const CartesianVector &nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    float initialGapSizeU(-10.f), maxGapSizeU(-10.f), maxProjectedGapSizeU(-10.f), hitLineDensityU(-10.f);

    LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pAlgorithm, nuVertexPosition, protoShowerU, TPC_VIEW_U, initialGapSizeU, maxGapSizeU, 
        maxProjectedGapSizeU, hitLineDensityU);

    float initialGapSizeV(-10.f), maxGapSizeV(-10.f), maxProjectedGapSizeV(-10.f), hitLineDensityV(-10.f);

    LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pAlgorithm, nuVertexPosition, protoShowerV, TPC_VIEW_V, initialGapSizeV, maxGapSizeV, 
        maxProjectedGapSizeV, hitLineDensityV);

    float initialGapSizeW(-10.f), maxGapSizeW(-10.f), maxProjectedGapSizeW(-10.f), hitLineDensityW(-10.f);

    LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pAlgorithm, nuVertexPosition, protoShowerW, TPC_VIEW_W, initialGapSizeW, maxGapSizeW, 
        maxProjectedGapSizeW, hitLineDensityW);

    float minInitialGapSize(-10.f), middleInitialGapSize(-10.f), maxInitialGapSize(-10.f);
    float minLargestGapSize(-10.f), middleLargestGapSize(-10.f), maxLargestGapSize(-10.f);
    float minLargestProjectedGapSize(-10.f), middleLargestProjectedGapSize(-10.f), maxLargestProjectedGapSize(-10.f);
    float minHitLineDensity(-10.f), middleHitLineDensity(-10.f), maxHitLineDensity(-10.f);

    LArConnectionPathwayHelper::GetMinMiddleMax(initialGapSizeU, initialGapSizeV, initialGapSizeW, minInitialGapSize, middleInitialGapSize, maxInitialGapSize);
    LArConnectionPathwayHelper::GetMinMiddleMax(maxGapSizeU, maxGapSizeV, maxGapSizeW, minLargestGapSize, middleLargestGapSize, maxLargestGapSize);
    LArConnectionPathwayHelper::GetMinMiddleMax(maxProjectedGapSizeU, maxProjectedGapSizeV, maxProjectedGapSizeW, minLargestProjectedGapSize, middleLargestProjectedGapSize, maxLargestProjectedGapSize);
    LArConnectionPathwayHelper::GetMinMiddleMax(hitLineDensityU, hitLineDensityV, hitLineDensityW, minHitLineDensity, middleHitLineDensity, maxHitLineDensity);

    electronTreeVariables.m_initialGapSizeW = initialGapSizeW;
    electronTreeVariables.m_minInitialGapSize = minInitialGapSize;
    electronTreeVariables.m_middleInitialGapSize = middleInitialGapSize;
    electronTreeVariables.m_maxInitialGapSize = std::min(maxInitialGapSize, 4.f);
    electronTreeVariables.m_largestGapSizeW = maxGapSizeW;
    electronTreeVariables.m_minLargestGapSize = minLargestGapSize;
    electronTreeVariables.m_middleLargestGapSize = middleLargestGapSize;
    electronTreeVariables.m_maxLargestGapSize = maxLargestGapSize;
    electronTreeVariables.m_largestProjectedGapSizeW = maxProjectedGapSizeW;
    electronTreeVariables.m_minLargestProjectedGapSize = std::min(minLargestProjectedGapSize, 2.f);
    electronTreeVariables.m_middleLargestProjectedGapSize = middleLargestProjectedGapSize;
    electronTreeVariables.m_maxLargestProjectedGapSize = maxLargestProjectedGapSize;
    electronTreeVariables.m_hitLineDensityW = hitLineDensityW;
    electronTreeVariables.m_minHitLineDensity = minHitLineDensity;
    electronTreeVariables.m_middleHitLineDensity = middleHitLineDensity;
    electronTreeVariables.m_maxHitLineDensity = maxHitLineDensity;

    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    try
    {
        // try 10 so i can make it consistent?
        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        CartesianVector startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        CartesianVector startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        CartesianVector startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
            protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
        {
                std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY - USING ANOTHER METHOD..." << std::endl;

            if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
                protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, protoShowerU.m_connectionPathway.m_startDirection, protoShowerV.m_connectionPathway.m_startDirection, 
                protoShowerW.m_connectionPathway.m_startDirection, connectionPathwayDirection3D))
            {
                std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
                return false;
            }
        }

        const CartesianVector projectedNuVertexU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_U));
        const CartesianVector projectedNuVertexV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_V));
        const CartesianVector projectedNuVertexW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_W));

        // The 3D hits of the connecting pathway might not exist -.- so we are going to have to try and make them...

        std::map<float, const CaloHit*> longitudinalProjectionsU, longitudinalProjectionsV, longitudinalProjectionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
            longitudinalProjectionsU[startDirectionU.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexU)] = pCaloHit;

        for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
            longitudinalProjectionsV[startDirectionV.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexV)] = pCaloHit;

        for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
            longitudinalProjectionsW[startDirectionW.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexW)] = pCaloHit;

        // maybe we want to limit this? to be like 10cm?
        const float length((middleShowerStart3D - nuVertexPosition).GetMagnitude());

        bool inPathway(true);
        float stepSize(2.0f);
        unsigned int count(0);

        // search for the closest point
        bool foundClosestPoint(false);
        CartesianVector closestPosition3D(0.f, 0.f, 0.f);

        ///////////////////////////
        //const CaloHit *pBestHitU(nullptr), *pBestHitV(nullptr), *pBestHitW(nullptr);
        ///////////////////////////

        // search for any gaps
        float gapSize(0.f);
        float distanceInGaps(0.f);
        float largestGapSize(-std::numeric_limits<float>::max());

        while (inPathway)
        {
            const CartesianVector lowerBoundary(nuVertexPosition + (connectionPathwayDirection3D * count * stepSize));
            const CartesianVector lowerBoundaryU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), lowerBoundary, TPC_VIEW_U));
            const CartesianVector lowerBoundaryV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), lowerBoundary, TPC_VIEW_V));
            const CartesianVector lowerBoundaryW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), lowerBoundary, TPC_VIEW_W));

            const float uLowerBoundaryProjectionL(startDirectionU.GetDotProduct(lowerBoundaryU - projectedNuVertexU));
            const float vLowerBoundaryProjectionL(startDirectionV.GetDotProduct(lowerBoundaryV - projectedNuVertexV));
            const float wLowerBoundaryProjectionL(startDirectionW.GetDotProduct(lowerBoundaryW - projectedNuVertexW));

            const CartesianVector upperBoundary(nuVertexPosition + (connectionPathwayDirection3D * ((count + 1) * stepSize)));
            const CartesianVector upperBoundaryU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), upperBoundary, TPC_VIEW_U));
            const CartesianVector upperBoundaryV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), upperBoundary, TPC_VIEW_V));
            const CartesianVector upperBoundaryW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), upperBoundary, TPC_VIEW_W));

            const float uUpperBoundaryProjectionL(startDirectionU.GetDotProduct(upperBoundaryU - projectedNuVertexU));
            const float vUpperBoundaryProjectionL(startDirectionV.GetDotProduct(upperBoundaryV - projectedNuVertexV));
            const float wUpperBoundaryProjectionL(startDirectionW.GetDotProduct(upperBoundaryW - projectedNuVertexW));

            //////////////////////////////
            /*
            for (const auto &entryU : longitudinalProjectionsU)
            {
                if ((entryU.first > uLowerBoundaryProjectionL) && (entryU.first < uUpperBoundaryProjectionL))
                {
                    const CartesianVector &hitPosition(entryU.second->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "hitPosition", GREEN, 2);
                }
            }

            for (const auto &entryV : longitudinalProjectionsV)
            {
                if ((entryV.first > vLowerBoundaryProjectionL) && (entryV.first < vUpperBoundaryProjectionL))
                {
                    const CartesianVector &hitPosition(entryV.second->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "hitPosition", GREEN, 2);
                }
            }

            for (const auto &entryW : longitudinalProjectionsW)
            {
                if ((entryW.first > wLowerBoundaryProjectionL) && (entryW.first < wUpperBoundaryProjectionL))
                {
                    const CartesianVector &hitPosition(entryW.second->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "hitPosition", GREEN, 2);
                }
            }

            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &lowerBoundaryU, "lowerBoundaryU", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &upperBoundaryU, "upperBoundaryU", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &lowerBoundaryV, "lowerBoundaryV", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &upperBoundaryV, "upperBoundaryV", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &lowerBoundaryW, "lowerBoundaryW", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &upperBoundaryW, "upperBoundaryW", BLACK, 2);
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            */
            //////////////////////////////

            bool foundHit3D(false);
            float lowestL(std::numeric_limits<float>::max());

            for (const auto &entryU : longitudinalProjectionsU)
            {
                if ((entryU.first > uLowerBoundaryProjectionL) && (entryU.first < uUpperBoundaryProjectionL))
                {
                    for (const auto &entryV : longitudinalProjectionsV)
                    {
                        if ((entryV.first > vLowerBoundaryProjectionL) && (entryV.first < vUpperBoundaryProjectionL))
                        {
                            for (const auto &entryW : longitudinalProjectionsW)
                            {
                                if ((entryW.first > wLowerBoundaryProjectionL) && (entryW.first < wUpperBoundaryProjectionL))
                                 {
                                     float metric(0.f);
                                     if (LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, entryU.second->GetPositionVector(), 
                                         entryV.second->GetPositionVector(), entryW.second->GetPositionVector(), 5.f, 2.f, metric))
                                     {
                                         foundHit3D = true;

                                         CartesianVector position3D(0.f, 0.f, 0.f);

                                         LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, entryU.second->GetPositionVector(), 
                                             entryV.second->GetPositionVector(), entryW.second->GetPositionVector(), position3D, metric);

                                         const float thisL(connectionPathwayDirection3D.GetDotProduct(position3D - nuVertexPosition));

                                         if (!foundClosestPoint && (thisL < lowestL))
                                         {
                                             lowestL = thisL;

                                             ///////////////////////////
                                             //pBestHitU = entryU.second; pBestHitV = entryV.second; pBestHitW = entryW.second;
                                             ///////////////////////////

                                             closestPosition3D = position3D;
                                         }
                                     }
                                 }
                            }
                        }
                    }
                }
            }

            if (!foundHit3D)
            {
                gapSize += stepSize;
            }
            else
            {
                if (gapSize > 1.f)
                    distanceInGaps += gapSize;

                if (gapSize > largestGapSize)
                    largestGapSize = gapSize;

                gapSize = 0.f;
            }

            if (!foundClosestPoint && foundHit3D)
            {
                ///////////////////////////
                /*
                const CartesianVector &closestHitU(pBestHitU->GetPositionVector());
                const CartesianVector &closestHitV(pBestHitV->GetPositionVector());
                const CartesianVector &closestHitW(pBestHitW->GetPositionVector());

                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestHitU, "closestHitU", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestHitV, "closestHitV", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestHitW, "closestHitW", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestPosition3D, "closestPosition3D", RED, 2);
                PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                */
                ///////////////////////////
                
                foundClosestPoint = true;
            }

            if (((upperBoundary - nuVertexPosition).GetMagnitude() > length) || (foundClosestPoint))
                inPathway = false;

            ++count;
        }

        if (foundClosestPoint)
        {
            const float distanceToNuVertex((closestPosition3D - nuVertexPosition).GetMagnitude());
            electronTreeVariables.m_initialRegionDistanceToNuVertex = distanceToNuVertex;
            //std::cout << "electronTreeVariables.m_initialRegionDistanceToNuVertex: " << electronTreeVariables.m_initialRegionDistanceToNuVertex << std::endl;
        }

        if (largestGapSize < 0.f)
            largestGapSize = -10.f;

        electronTreeVariables.m_initialRegionDistanceInGaps = distanceInGaps;
        //std::cout << "electronTreeVariables.m_initialRegionDistanceInGaps: " << electronTreeVariables.m_initialRegionDistanceInGaps << std::endl;
        electronTreeVariables.m_initialRegionMaxGapSize = largestGapSize;
        //std::cout << "electronTreeVariables.m_initialRegionMaxGapSize: " << electronTreeVariables.m_initialRegionMaxGapSize << std::endl;
    }
    catch (...)
    {
        //std::cout << "am i here??" << std::endl;
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// I think I can do this at the hit level and avoid outliers because i'm already following a fit?    
bool LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pandora::Algorithm *const pAlgorithm, const CartesianVector &nuVertexPosition, const ProtoShower &protoShower, 
    const HitType hitType, float &initialGapSize, float &maxGapSize, float &maxProjectedGapSize, float &hitLineDensity)
{
    // First look at longitudinal projections
    ///////////////////////////////////////////////////////

    maxProjectedGapSize = -10.f;

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, hitType));
    const CartesianVector &startDirection(protoShower.m_connectionPathway.m_startDirection);

    FloatVector longitudinalProjections;

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
        longitudinalProjections.push_back(startDirection.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexPosition));

    std::sort(longitudinalProjections.begin(), longitudinalProjections.end());

    const long unsigned int nSampleHits(10);
    const unsigned int nIterations(std::min(longitudinalProjections.size(), nSampleHits) - 1); // basically how many hits you want to consider

    for (unsigned int i = 0; i < nIterations; ++i)
        maxProjectedGapSize = std::max(std::fabs(longitudinalProjections[i] - longitudinalProjections[i + 1]), maxProjectedGapSize);

    // average gap between hits
    hitLineDensity = std::fabs(longitudinalProjections[nIterations] - longitudinalProjections.front()) / static_cast<float>(nIterations);

    // Next look at absolute disance between hits
    ///////////////////////////////////////////////////////

    CaloHitVector spineHitVector(protoShower.m_spineHitList.begin(), protoShower.m_spineHitList.end());
    std::sort(spineHitVector.begin(), spineHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(projectedNuVertexPosition));

    initialGapSize = (projectedNuVertexPosition - spineHitVector.front()->GetPositionVector()).GetMagnitude();

    spineHitVector.resize(nIterations);

    maxGapSize = -10.f;

    for (unsigned int i = 0; i < spineHitVector.size(); ++i)
    {
        float closestDistance(std::numeric_limits<float>::max());

        for (unsigned int j = 0; j < spineHitVector.size(); ++j)
        {
            if (i == j)
                continue;

            float distance((spineHitVector[i]->GetPositionVector() - spineHitVector[j]->GetPositionVector()).GetMagnitude());

            if (distance < closestDistance)
                closestDistance = distance;
        }

        if (closestDistance > maxGapSize)
            maxGapSize = closestDistance;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillPostShowerStartVariables(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const CartesianVector &nuVertexPosition, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    try
    {
        // try 10 so i can make it consistent?
        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        ///////////////////////////////////////////////////
        const LayerFitResultMap &layerFitResultMapU(spineFitU.GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapV(spineFitV.GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapW(spineFitW.GetLayerFitResultMap());

        const CartesianVector &startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection());
        const CartesianVector &startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection());
        const CartesianVector &startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection());

        const CartesianVector projectedShowerStartU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_U));
        const CartesianVector projectedShowerStartV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_V));
        const CartesianVector projectedShowerStartW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_W));

        const CartesianVector projectedNuVertexU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_U));
        const CartesianVector projectedNuVertexV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_V));
        const CartesianVector projectedNuVertexW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_W));

        // attempt to get a 3D direction...
        float lShowerStartU(0.f), lShowerStartV(0.f), lShowerStartW(0.f);
        float tShowerStartU(0.f), tShowerStartV(0.f), tShowerStartW(0.f);

        spineFitU.GetLocalPosition(projectedShowerStartU, lShowerStartU, tShowerStartU);
        spineFitV.GetLocalPosition(projectedShowerStartV, lShowerStartV, tShowerStartV);
        spineFitW.GetLocalPosition(projectedShowerStartW, lShowerStartW, tShowerStartW);

        const int showerStartLayerU(spineFitU.GetLayer(lShowerStartU));
        int closestLayerU(std::numeric_limits<int>::max());

        for (const auto &entry : layerFitResultMapU)
        {
            if (std::fabs(entry.first - showerStartLayerU) < std::fabs(entry.first - closestLayerU))
                closestLayerU = entry.first;
        } 

        const int showerStartLayerV(spineFitV.GetLayer(lShowerStartV));
        int closestLayerV(std::numeric_limits<int>::max());

        for (const auto &entry : layerFitResultMapV)
        {
            if (std::fabs(entry.first - showerStartLayerV) < std::fabs(entry.first - closestLayerV))
                closestLayerV = entry.first;
        } 

        const int showerStartLayerW(spineFitW.GetLayer(lShowerStartW));
        int closestLayerW(std::numeric_limits<int>::max());

        for (const auto &entry : layerFitResultMapW)
        {
            if (std::fabs(entry.first - showerStartLayerW) < std::fabs(entry.first - closestLayerW))
                closestLayerW = entry.first;
        } 

        if ((layerFitResultMapU.find(closestLayerU) == layerFitResultMapU.end()) || (layerFitResultMapV.find(closestLayerV) == layerFitResultMapV.end()) ||
            (layerFitResultMapW.find(closestLayerW) == layerFitResultMapW.end()))
        {
            ////std::cout << "shower start layer is not in layer fit result map" << std::endl;
            return false;
        }

        const float gradientU(layerFitResultMapU.at(closestLayerU).GetGradient());
        const float gradientV(layerFitResultMapV.at(closestLayerV).GetGradient());
        const float gradientW(layerFitResultMapW.at(closestLayerW).GetGradient());
        CartesianVector showerDirectionU(0.f, 0.f, 0.f), showerDirectionV(0.f, 0.f, 0.f), showerDirectionW(0.f, 0.f, 0.f);

        spineFitU.GetGlobalDirection(gradientU, showerDirectionU);
        spineFitV.GetGlobalDirection(gradientV, showerDirectionV);
        spineFitW.GetGlobalDirection(gradientW, showerDirectionW);

        int nHitsU(-10);
        float scatterAngleU(-10.f), openingAngleU(-10.f), openingAngleAsymmetryU(-10.f), nuVertexHitAsymmetryU(-10.f), nuVertexEnergyAsymmetryU(-0.5f),
            showerStartHitAsymmetryU(-10.f), showerStartEnergyAsymmetryU(-0.5f), nuVertexMeanRadialDistanceU(-10.f), nuVertexEnergyWeightedMeanRadialDistanceU(-10.f),
            showerStartMeanRadialDistanceU(-10.f), showerStartEnergyWeightedMeanRadialDistanceU(-10.f), nuVertexMoliereRadiusU(-10.f), showerStartMoliereRadiusU(-10.f), 
            positiveOpeningAngleU(-10.f), negativeOpeningAngleU(-10.f), showerApexLU(-9999.f), showerApexTU(-10.f), fitShowerStartLU(-10.f), fitShowerStartTU(-10.f), 
            foundHitRatioU(-0.5f);

        LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pAlgorithm, protoShowerU.m_spineHitList, pShowerPfo, TPC_VIEW_U, isDownstream, projectedNuVertexU, 
            projectedShowerStartU, showerDirectionU, startDirectionU, nHitsU, scatterAngleU, openingAngleU, openingAngleAsymmetryU, nuVertexHitAsymmetryU, nuVertexEnergyAsymmetryU, 
            showerStartHitAsymmetryU, showerStartEnergyAsymmetryU, nuVertexMeanRadialDistanceU, nuVertexEnergyWeightedMeanRadialDistanceU, showerStartMeanRadialDistanceU, 
            showerStartEnergyWeightedMeanRadialDistanceU, nuVertexMoliereRadiusU, showerStartMoliereRadiusU, positiveOpeningAngleU, negativeOpeningAngleU, 
            showerApexLU, showerApexTU, fitShowerStartLU, fitShowerStartTU, foundHitRatioU);

        openingAngleU = std::max(positiveOpeningAngleU, negativeOpeningAngleU);

        int nHitsV(-10);
        float scatterAngleV(-10.f), openingAngleV(-10.f), openingAngleAsymmetryV(-10.f), nuVertexHitAsymmetryV(-10.f), nuVertexEnergyAsymmetryV(-0.5f),
            showerStartHitAsymmetryV(-10.f), showerStartEnergyAsymmetryV(-0.5f), nuVertexMeanRadialDistanceV(-10.f), nuVertexEnergyWeightedMeanRadialDistanceV(-10.f),
            showerStartMeanRadialDistanceV(-10.f), showerStartEnergyWeightedMeanRadialDistanceV(-10.f), nuVertexMoliereRadiusV(-10.f), showerStartMoliereRadiusV(-10.f),
            positiveOpeningAngleV(-10.f), negativeOpeningAngleV(-10.f), showerApexLV(-9999.f), showerApexTV(-10.f), fitShowerStartLV(-10.f), fitShowerStartTV(-10.f),
            foundHitRatioV(-0.5f);

        LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pAlgorithm, protoShowerV.m_spineHitList, pShowerPfo, TPC_VIEW_V, isDownstream, projectedNuVertexV, 
            projectedShowerStartV, showerDirectionV, startDirectionV, nHitsV, scatterAngleV, openingAngleV, openingAngleAsymmetryV, nuVertexHitAsymmetryV, nuVertexEnergyAsymmetryV, 
            showerStartHitAsymmetryV, showerStartEnergyAsymmetryV, nuVertexMeanRadialDistanceV, nuVertexEnergyWeightedMeanRadialDistanceV, showerStartMeanRadialDistanceV, 
            showerStartEnergyWeightedMeanRadialDistanceV, nuVertexMoliereRadiusV, showerStartMoliereRadiusV, positiveOpeningAngleV, negativeOpeningAngleV, 
            showerApexLV, showerApexTV, fitShowerStartLV, fitShowerStartTV, foundHitRatioV);

        openingAngleV = std::max(positiveOpeningAngleV, negativeOpeningAngleV);

        int nHitsW(-10);
        float scatterAngleW(-10.f), openingAngleW(-10.f), openingAngleAsymmetryW(-10.f), nuVertexHitAsymmetryW(-10.f), nuVertexEnergyAsymmetryW(-0.5f),
            showerStartHitAsymmetryW(-10.f), showerStartEnergyAsymmetryW(-0.5f), nuVertexMeanRadialDistanceW(-10.f), nuVertexEnergyWeightedMeanRadialDistanceW(-10.f),
            showerStartMeanRadialDistanceW(-10.f), showerStartEnergyWeightedMeanRadialDistanceW(-10.f), nuVertexMoliereRadiusW(-10.f), showerStartMoliereRadiusW(-10.f),
            positiveOpeningAngleW(-10.f), negativeOpeningAngleW(-10.f), showerApexLW(-9999.f), showerApexTW(-10.f), fitShowerStartLW(-10.f), fitShowerStartTW(-10.f),
            foundHitRatioW(-0.5f);

        LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pAlgorithm, protoShowerW.m_spineHitList, pShowerPfo, TPC_VIEW_W, isDownstream, projectedNuVertexW, 
            projectedShowerStartW, showerDirectionW, startDirectionW, nHitsW, scatterAngleW, openingAngleW, openingAngleAsymmetryW, nuVertexHitAsymmetryW, nuVertexEnergyAsymmetryW, 
            showerStartHitAsymmetryW, showerStartEnergyAsymmetryW, nuVertexMeanRadialDistanceW, nuVertexEnergyWeightedMeanRadialDistanceW, showerStartMeanRadialDistanceW, 
            showerStartEnergyWeightedMeanRadialDistanceW, nuVertexMoliereRadiusW, showerStartMoliereRadiusW, positiveOpeningAngleW, negativeOpeningAngleW, 
            showerApexLW, showerApexTW, fitShowerStartLW, fitShowerStartTW, foundHitRatioW);

        openingAngleW = std::max(positiveOpeningAngleW, negativeOpeningAngleW);

        float minNHits(-10.f), middleNHits(-10.f), maxNHits(-10.f);
        float minScatterAngle(-10.f), middleScatterAngle(-10.f), maxScatterAngle(-10.f);
        float minOpeningAngle(-10.f), middleOpeningAngle(-10.f), maxOpeningAngle(-10.f);
        float minOpeningAngleAsymmetry(-10.f), middleOpeningAngleAsymmetry(-10.f), maxOpeningAngleAsymmetry(-10.f);
        float minNuVertexHitAsymmetry(-10.f), middleNuVertexHitAsymmetry(-10.f), maxNuVertexHitAsymmetry(-10.f);
        float minNuVertexEnergyAsymmetry(-0.5f), middleNuVertexEnergyAsymmetry(-0.5f), maxNuVertexEnergyAsymmetry(-0.5f);
        float minShowerStartHitAsymmetry(-10.f), middleShowerStartHitAsymmetry(-10.f), maxShowerStartHitAsymmetry(-10.f);
        float minShowerStartEnergyAsymmetry(-0.5f), middleShowerStartEnergyAsymmetry(-0.5f), maxShowerStartEnergyAsymmetry(-0.5f);
        float minNuVertexMeanRadialDistance(-10.f), middleNuVertexMeanRadialDistance(-10.f), maxNuVertexMeanRadialDistance(-10.f);
        float minNuVertexEnergyWeightedMeanRadialDistance(-10.f), middleNuVertexEnergyWeightedMeanRadialDistance(-10.f), maxNuVertexEnergyWeightedMeanRadialDistance(-10.f);
        float minShowerStartMeanRadialDistance(-10.f), middleShowerStartMeanRadialDistance(-10.f), maxShowerStartMeanRadialDistance(-10.f);
        float minShowerStartEnergyWeightedMeanRadialDistance(-10.f), middleShowerStartEnergyWeightedMeanRadialDistance(-10.f), maxShowerStartEnergyWeightedMeanRadialDistance(-10.f);
        float minNuVertexMoliereRadius(-10.f), middleNuVertexMoliereRadius(-10.f), maxNuVertexMoliereRadius(-10.f);
        float minShowerStartMoliereRadius(-10.f), middleShowerStartMoliereRadius(-10.f), maxShowerStartMoliereRadius(-10.f);
        float minShowerApexL(-9999.f), middleShowerApexL(-9999.f), maxShowerApexL(-9999.f);
        float minShowerApexT(-10.f), middleShowerApexT(-10.f), maxShowerApexT(-10.f);
        float minFoundHitRatio(-0.5f), middleFoundHitRatio(-0.5), maxFoundHitRatio(-0.5f);

        LArConnectionPathwayHelper::GetMinMiddleMax(nHitsU, nHitsV, nHitsW, minNHits, middleNHits, maxNHits);
        LArConnectionPathwayHelper::GetMinMiddleMax(scatterAngleU, scatterAngleV, scatterAngleW, minScatterAngle, middleScatterAngle, maxScatterAngle);
        LArConnectionPathwayHelper::GetMinMiddleMax(openingAngleU, openingAngleV, openingAngleW, minOpeningAngle, middleOpeningAngle, maxOpeningAngle);
        LArConnectionPathwayHelper::GetMinMiddleMax(openingAngleAsymmetryU, openingAngleAsymmetryV, openingAngleAsymmetryW, minOpeningAngleAsymmetry, 
            middleOpeningAngleAsymmetry, maxOpeningAngleAsymmetry);
        LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexHitAsymmetryU, nuVertexHitAsymmetryV, nuVertexHitAsymmetryW, minNuVertexHitAsymmetry, 
            middleNuVertexHitAsymmetry, maxNuVertexHitAsymmetry);
        LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexEnergyAsymmetryU, nuVertexEnergyAsymmetryV, nuVertexEnergyAsymmetryW, minNuVertexEnergyAsymmetry, 
            middleNuVertexEnergyAsymmetry, maxNuVertexEnergyAsymmetry);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerStartHitAsymmetryU, showerStartHitAsymmetryV, showerStartHitAsymmetryW, minShowerStartHitAsymmetry, 
            middleShowerStartHitAsymmetry, maxShowerStartHitAsymmetry);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerStartEnergyAsymmetryU, showerStartEnergyAsymmetryV, showerStartEnergyAsymmetryW, minShowerStartEnergyAsymmetry, 
            middleShowerStartEnergyAsymmetry, maxShowerStartEnergyAsymmetry);
        LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexMeanRadialDistanceU, nuVertexMeanRadialDistanceV, nuVertexMeanRadialDistanceW, minNuVertexMeanRadialDistance, 
            middleNuVertexMeanRadialDistance, maxNuVertexMeanRadialDistance);
        LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexEnergyWeightedMeanRadialDistanceU, nuVertexEnergyWeightedMeanRadialDistanceV, nuVertexEnergyWeightedMeanRadialDistanceW, 
            minNuVertexEnergyWeightedMeanRadialDistance, middleNuVertexEnergyWeightedMeanRadialDistance, maxNuVertexEnergyWeightedMeanRadialDistance);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerStartMeanRadialDistanceU, showerStartMeanRadialDistanceV, showerStartMeanRadialDistanceW, minShowerStartMeanRadialDistance, 
            middleShowerStartMeanRadialDistance, maxShowerStartMeanRadialDistance);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerStartEnergyWeightedMeanRadialDistanceU, showerStartEnergyWeightedMeanRadialDistanceV, showerStartEnergyWeightedMeanRadialDistanceW, 
            minShowerStartEnergyWeightedMeanRadialDistance, middleShowerStartEnergyWeightedMeanRadialDistance, maxShowerStartEnergyWeightedMeanRadialDistance);
        LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexMoliereRadiusU, nuVertexMoliereRadiusV, nuVertexMoliereRadiusW, minNuVertexMoliereRadius, 
            middleNuVertexMoliereRadius, maxNuVertexMoliereRadius);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerStartMoliereRadiusU, showerStartMoliereRadiusV, showerStartMoliereRadiusW, minShowerStartMoliereRadius, 
            middleShowerStartMoliereRadius, maxShowerStartMoliereRadius);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerApexLU, showerApexLV, showerApexLW, minShowerApexL, middleShowerApexL, maxShowerApexL);
        LArConnectionPathwayHelper::GetMinMiddleMax(showerApexTU, showerApexTV, showerApexTW, minShowerApexT, middleShowerApexT, maxShowerApexT);
        LArConnectionPathwayHelper::GetMinMiddleMax(foundHitRatioU, foundHitRatioV, foundHitRatioW, minFoundHitRatio, middleFoundHitRatio, maxFoundHitRatio);

        electronTreeVariables.m_postShowerStartNHitsW = nHitsW;
        electronTreeVariables.m_minNPostShowerStartHits = minNHits;
        electronTreeVariables.m_middleNPostShowerStartHits = middleNHits;
        electronTreeVariables.m_maxNPostShowerStartHits = std::min(maxNHits, 2000.f);
        electronTreeVariables.m_postShowerStartScatterAngleW = scatterAngleW;
        electronTreeVariables.m_minPostShowerStartScatterAngle = minScatterAngle;
        electronTreeVariables.m_middlePostShowerStartScatterAngle = middleScatterAngle;
        electronTreeVariables.m_maxPostShowerStartScatterAngle = std::min(maxScatterAngle, 40.f);
        electronTreeVariables.m_postShowerStartOpeningAngleW = std::min(openingAngleW, 20.f);
        electronTreeVariables.m_minPostShowerStartOpeningAngle = std::min(minOpeningAngle, 20.f);
        electronTreeVariables.m_middlePostShowerStartOpeningAngle = std::min(middleOpeningAngle, 20.f);
        electronTreeVariables.m_maxPostShowerStartOpeningAngle = std::min(maxOpeningAngle, 20.f);
        electronTreeVariables.m_postShowerStartOpeningAngleAsymmetryW = openingAngleAsymmetryW;
        electronTreeVariables.m_minPostShowerStartOpeningAngleAsymmetry = minOpeningAngleAsymmetry;
        electronTreeVariables.m_middlePostShowerStartOpeningAngleAsymmetry = middleOpeningAngleAsymmetry;
        electronTreeVariables.m_maxPostShowerStartOpeningAngleAsymmetry = maxOpeningAngleAsymmetry;
        electronTreeVariables.m_postShowerStartNuVertexHitAsymmetryW = nuVertexHitAsymmetryW;
        electronTreeVariables.m_minPostShowerStartNuVertexHitAsymmetry = minNuVertexHitAsymmetry;
        electronTreeVariables.m_middlePostShowerStartNuVertexHitAsymmetry = middleNuVertexHitAsymmetry;
        electronTreeVariables.m_maxPostShowerStartNuVertexHitAsymmetry = maxNuVertexHitAsymmetry;
        electronTreeVariables.m_postShowerStartNuVertexEnergyAsymmetryW = nuVertexEnergyAsymmetryW;
        electronTreeVariables.m_minPostShowerStartNuVertexEnergyAsymmetry = minNuVertexEnergyAsymmetry;
        electronTreeVariables.m_middlePostShowerStartNuVertexEnergyAsymmetry = middleNuVertexEnergyAsymmetry;
        electronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry = maxNuVertexEnergyAsymmetry;
        electronTreeVariables.m_postShowerStartShowerStartHitAsymmetryW = showerStartHitAsymmetryW;
        electronTreeVariables.m_minPostShowerStartShowerStartHitAsymmetry = minShowerStartHitAsymmetry;
        electronTreeVariables.m_middlePostShowerStartShowerStartHitAsymmetry = middleShowerStartHitAsymmetry;
        electronTreeVariables.m_maxPostShowerStartShowerStartHitAsymmetry = maxShowerStartHitAsymmetry;
        electronTreeVariables.m_postShowerStartShowerStartEnergyAsymmetryW = showerStartEnergyAsymmetryW;
        electronTreeVariables.m_minPostShowerStartShowerStartEnergyAsymmetry = minShowerStartEnergyAsymmetry;
        electronTreeVariables.m_middlePostShowerStartShowerStartEnergyAsymmetry = middleShowerStartEnergyAsymmetry;
        electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry = maxShowerStartEnergyAsymmetry;
        electronTreeVariables.m_postShowerStartNuVertexMeanRadialDistanceW = nuVertexMeanRadialDistanceW;
        electronTreeVariables.m_minPostShowerStartNuVertexMeanRadialDistance = minNuVertexMeanRadialDistance;
        electronTreeVariables.m_middlePostShowerStartNuVertexMeanRadialDistance = middleNuVertexMeanRadialDistance;
        electronTreeVariables.m_maxPostShowerStartNuVertexMeanRadialDistance = maxNuVertexMeanRadialDistance;
        electronTreeVariables.m_postShowerStartNuVertexEnergyWeightedMeanRadialDistanceW = nuVertexEnergyWeightedMeanRadialDistanceW;
        electronTreeVariables.m_minPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = minNuVertexEnergyWeightedMeanRadialDistance;
        electronTreeVariables.m_middlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance = middleNuVertexEnergyWeightedMeanRadialDistance;
        electronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = std::min(maxNuVertexEnergyWeightedMeanRadialDistance, 20.f);
        electronTreeVariables.m_postShowerStartShowerStartMeanRadialDistanceW = showerStartMeanRadialDistanceW;
        electronTreeVariables.m_minPostShowerStartShowerStartMeanRadialDistance = minShowerStartMeanRadialDistance;
        electronTreeVariables.m_middlePostShowerStartShowerStartMeanRadialDistance = middleShowerStartMeanRadialDistance;
        electronTreeVariables.m_maxPostShowerStartShowerStartMeanRadialDistance = maxShowerStartMeanRadialDistance;
        electronTreeVariables.m_postShowerStartShowerStartEnergyWeightedMeanRadialDistanceW = showerStartEnergyWeightedMeanRadialDistanceW;
        electronTreeVariables.m_minPostShowerStartShowerStartEnergyWeightedMeanRadialDistance = minShowerStartEnergyWeightedMeanRadialDistance;
        electronTreeVariables.m_middlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance = middleShowerStartEnergyWeightedMeanRadialDistance;
        electronTreeVariables.m_maxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance = maxShowerStartEnergyWeightedMeanRadialDistance;
        electronTreeVariables.m_postShowerStartNuVertexMoliereRadiusW = nuVertexMoliereRadiusW;
        electronTreeVariables.m_minPostShowerStartNuVertexMoliereRadius = minNuVertexMoliereRadius;
        electronTreeVariables.m_middlePostShowerStartNuVertexMoliereRadius = middleNuVertexMoliereRadius;
        electronTreeVariables.m_maxPostShowerStartNuVertexMoliereRadius = maxNuVertexMoliereRadius;
        electronTreeVariables.m_postShowerStartShowerStartMoliereRadiusW = showerStartMoliereRadiusW;
        electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius = std::min(minShowerStartMoliereRadius, 10.f);
        electronTreeVariables.m_middlePostShowerStartShowerStartMoliereRadius = middleShowerStartMoliereRadius;
        electronTreeVariables.m_maxPostShowerStartShowerStartMoliereRadius = maxShowerStartMoliereRadius;
        electronTreeVariables.m_positiveOpeningAngleW = positiveOpeningAngleW;
        electronTreeVariables.m_negativeOpeningAngleW = negativeOpeningAngleW;
        //electronTreeVariables.m_maxOpeningAngleW = std::min(std::max(positiveOpeningAngleW, negativeOpeningAngleW), 20.f);
        electronTreeVariables.m_showerApexLW = showerApexLW;
        electronTreeVariables.m_minShowerApexL = minShowerApexL;
        electronTreeVariables.m_middleShowerApexL = middleShowerApexL;
        electronTreeVariables.m_maxShowerApexL = maxShowerApexL;
        electronTreeVariables.m_showerApexTW = showerApexTW;
        electronTreeVariables.m_minShowerApexT = minShowerApexT;
        electronTreeVariables.m_middleShowerApexT = middleShowerApexT;
        electronTreeVariables.m_maxShowerApexT = maxShowerApexT;
        electronTreeVariables.m_foundHitRatioW = foundHitRatioW;
        electronTreeVariables.m_minFoundHitRatio = minFoundHitRatio;
        electronTreeVariables.m_middleFoundHitRatio = middleFoundHitRatio;
        electronTreeVariables.m_maxFoundHitRatio = std::min(maxFoundHitRatio, 1.5f);
        electronTreeVariables.m_fitShowerStartLW = fitShowerStartLW;
        electronTreeVariables.m_fitShowerStartTW = fitShowerStartTW;

        // Now do 3D shower variables....
        //////////////////////////////////////////////////

        CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
            protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
        {
            if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
                protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, protoShowerU.m_connectionPathway.m_startDirection, protoShowerV.m_connectionPathway.m_startDirection, 
                protoShowerW.m_connectionPathway.m_startDirection, connectionPathwayDirection3D))
            {
                //std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
                return false;
            }
        }

        ////////////////////////////////////
        /*
        const CartesianVector endU(projectedShowerStartU + (showerDirectionU * 100.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedShowerStartU, &endU, "PROJECTION U", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "U SPINE1", GREEN);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        const CartesianVector endV(projectedShowerStartV + (showerDirectionV * 100.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedShowerStartV, &endV, "PROJECTION V", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "V SPINE1", GREEN);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        const CartesianVector endW(projectedShowerStartW + (showerDirectionW * 100.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedShowerStartW, &endW, "PROJECTION W", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "W SPINE1", GREEN);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        ////////////////////////////////////

        CartesianVector postShowerStartDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, projectedShowerStartU, projectedShowerStartV, projectedShowerStartW, middleShowerStart3D,
            showerDirectionU, showerDirectionV, showerDirectionW, postShowerStartDirection3D))
        {
            //std::cout << "CANNOT FIND 3D POST SHOWER START DIRECTION... SO I AM USING THE CONNECTION PATHWAY..." << std::endl;
            postShowerStartDirection3D = connectionPathwayDirection3D;
        }

        CaloHitList caloHitList3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHitList3D);

        CaloHitList postShowerHitList3D; // Do i use all hits?? or all hits in event?? o.O 
        CartesianPointVector postShowerPositions3D;

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - middleShowerStart3D);
            const float l(postShowerStartDirection3D.GetDotProduct(displacement));

            const float openingAngle(displacement.GetOpeningAngle(postShowerStartDirection3D) * 180.f / M_PI);

            if ((l > 0.f) && (openingAngle < 15.f))
            {
                postShowerHitList3D.push_back(pCaloHit);
                postShowerPositions3D.push_back(pCaloHit->GetPositionVector());
            }
        }

        int nPostShowerHits(postShowerHitList3D.size());
        electronTreeVariables.m_postShowerStartNHits = nPostShowerHits;
        //std::cout << "electronTreeVariables.m_postShowerStartNHits: " << electronTreeVariables.m_postShowerStartNHits << std::endl;

        ThreeDSlidingFitResult slidingFitResult3D(&postShowerPositions3D, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            if (std::find(postShowerHitList3D.begin(), postShowerHitList3D.end(), pCaloHit) != postShowerHitList3D.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - middleShowerStart3D);
            const float l(slidingFitResult3D.GetGlobalMinLayerDirection().GetDotProduct(displacement));

            float openingAngle(displacement.GetOpeningAngle(slidingFitResult3D.GetGlobalMinLayerDirection()) * 180.f / M_PI);
            openingAngle = isDownstream ? openingAngle : 360.f - openingAngle;

            if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (openingAngle < 25.f))
            {
                postShowerHitList3D.push_back(pCaloHit);
                postShowerPositions3D.push_back(pCaloHit->GetPositionVector());
            }
        }

        ///////////////////////////        
        /*
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &postShowerHitList3D, "postShowerHitList3D", RED);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        ///////////////////////////        

        ThreeDSlidingFitResult fullSlidingFitResult3D(&postShowerPositions3D, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const float length((fullSlidingFitResult3D.GetGlobalMaxLayerPosition() - fullSlidingFitResult3D.GetGlobalMinLayerPosition()).GetMagnitude());
        electronTreeVariables.m_postShowerStartLength = length;
        //std::cout << "electronTreeVariables.m_postShowerStartLength: " << electronTreeVariables.m_postShowerStartLength << std::endl;

        nPostShowerHits = postShowerHitList3D.size();
        electronTreeVariables.m_postShowerStartNHits = nPostShowerHits;
        //std::cout << "electronTreeVariables.m_postShowerStartNHits (after extra hit addition): " << electronTreeVariables.m_postShowerStartNHits << std::endl;

        float pathwayShowerAngle(fullSlidingFitResult3D.GetGlobalMinLayerDirection().GetOpeningAngle(connectionPathwayDirection3D) * 180.f / M_PI);
        pathwayShowerAngle = isDownstream ? pathwayShowerAngle : (360.f - pathwayShowerAngle);
        electronTreeVariables.m_postShowerStartScatterAngle = pathwayShowerAngle;
        //std::cout << "electronTreeVariables.m_postShowerStartScatterAngle: " << electronTreeVariables.m_postShowerStartScatterAngle << std::endl; 
 
        try
        {
            const CartesianVector edgeDirectionAxis(fullSlidingFitResult3D.GetGlobalMinLayerDirection());
            CartesianVector edgeOrthoAxis1(edgeDirectionAxis.GetCrossProduct(CartesianVector(1.f, 1.f, 1.f)).GetUnitVector());
            const CartesianVector edgeOrthoAxis2(edgeDirectionAxis.GetCrossProduct(edgeOrthoAxis1).GetUnitVector());

            std::map<int, float> positiveEdges, negativeEdges;

            for (const CaloHit *const pCaloHit : postShowerHitList3D)
            {
                const CartesianVector position(pCaloHit->GetPositionVector() - middleShowerStart3D);

                const float thisT(edgeDirectionAxis.GetCrossProduct(position).GetMagnitude());
                const float thisL(edgeDirectionAxis.GetDotProduct(position));
                const float orthoL(edgeOrthoAxis1.GetDotProduct(position));

                std::map<int, float> &edgeMap(orthoL > 0.f ? positiveEdges : negativeEdges);

                const int lIndex(std::floor(thisL / 2.f));

                edgeMap[lIndex] = (edgeMap.find(lIndex) == edgeMap.end() ? thisT : std::max(edgeMap[lIndex] , thisT));;
            }

            CartesianPointVector positiveEdgePositions, negativeEdgePositions;

            for (auto &entry : positiveEdges)
                positiveEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

            for (auto &entry : negativeEdges)
                negativeEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

            const TwoDSlidingFitResult positiveEdgeFit(&positiveEdgePositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
            const TwoDSlidingFitResult negativeEdgeFit(&negativeEdgePositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

            const CartesianVector positiveMinLayer(positiveEdgeFit.GetGlobalMinLayerPosition());
            const CartesianVector positiveMaxLayer(positiveEdgeFit.GetGlobalMaxLayerPosition());
            const CartesianVector negativeMinLayer(negativeEdgeFit.GetGlobalMinLayerPosition());
            const CartesianVector negativeMaxLayer(negativeEdgeFit.GetGlobalMaxLayerPosition());

            const float positiveGradient((positiveMaxLayer.GetZ() - positiveMinLayer.GetZ()) / (positiveMaxLayer.GetX() - positiveMinLayer.GetX()));
            const float negativeGradient((negativeMaxLayer.GetZ() - negativeMinLayer.GetZ()) / (negativeMaxLayer.GetX() - negativeMinLayer.GetX()));
            const float positiveIntercept(positiveMaxLayer.GetZ() - (positiveGradient * positiveMaxLayer.GetX()));
            const float negativeIntercept(negativeMaxLayer.GetZ() - (negativeGradient * negativeMaxLayer.GetX()));

            const float collisionT((-1.f) * (positiveIntercept + negativeIntercept) / (positiveGradient + negativeGradient));
            const float collisionL((positiveGradient * collisionT) + positiveIntercept);

            const float positiveStartL(positiveIntercept);
            const float negativeStartL(negativeIntercept);

            const CartesianVector globalPositiveMinLayer(middleShowerStart3D + (edgeDirectionAxis * positiveMinLayer.GetZ()) + (edgeOrthoAxis1 * positiveMinLayer.GetX()));
            const CartesianVector globalPositiveMaxLayer(middleShowerStart3D + (edgeDirectionAxis * positiveMaxLayer.GetZ()) + (edgeOrthoAxis1 * positiveMaxLayer.GetX()));
            const CartesianVector globalNegativeMinLayer(middleShowerStart3D + (edgeDirectionAxis * negativeMinLayer.GetZ()) + (edgeOrthoAxis1 * negativeMinLayer.GetX()));
            const CartesianVector globalNegativeMaxLayer(middleShowerStart3D + (edgeDirectionAxis * negativeMaxLayer.GetZ()) + (edgeOrthoAxis1 * negativeMaxLayer.GetX()));
            const CartesianVector collisionPoint(middleShowerStart3D + (edgeDirectionAxis * collisionL));
            const CartesianVector positiveStart(middleShowerStart3D + (edgeDirectionAxis * positiveStartL));
            const CartesianVector negativeStart(middleShowerStart3D + (edgeDirectionAxis * negativeStartL));

            const CartesianVector positiveEdgeVector((globalPositiveMaxLayer - globalPositiveMinLayer).GetUnitVector());
            const CartesianVector negativeEdgeVector((globalNegativeMaxLayer - globalNegativeMinLayer).GetUnitVector());

            const float positiveOpeningAngle(edgeDirectionAxis.GetOpeningAngle(positiveEdgeVector) * 180.f / M_PI);
            const float negativeOpeningAngle(edgeDirectionAxis.GetOpeningAngle(negativeEdgeVector) * 180.f / M_PI);

            const float minOpeningAngle3D(std::min(positiveOpeningAngle, negativeOpeningAngle));
            electronTreeVariables.m_postShowerStartMinHalfOpeningAngle = minOpeningAngle3D;
            //std::cout << "electronTreeVariables.m_postShowerStartMinHalfOpeningAngle: " << electronTreeVariables.m_postShowerStartMinHalfOpeningAngle << std::endl;

            const float maxOpeningAngle3D(std::max(positiveOpeningAngle, negativeOpeningAngle));
            electronTreeVariables.m_postShowerStartMaxHalfOpeningAngle = maxOpeningAngle3D;
            //std::cout << "electronTreeVariables.m_postShowerStartMaxHalfOpeningAngle: " << electronTreeVariables.m_postShowerStartMaxHalfOpeningAngle << std::endl;

            const float openingAngle3D(positiveOpeningAngle + negativeOpeningAngle);
            electronTreeVariables.m_postShowerStartOpeningAngle = openingAngle3D;
            //std::cout << "electronTreeVariables.m_postShowerStartOpeningAngle: " << electronTreeVariables.m_postShowerStartOpeningAngle << std::endl;

            const float openingAngleAsymmetry3D(std::fabs(positiveOpeningAngle - negativeOpeningAngle));
            electronTreeVariables.m_postShowerStartOpeningAngleAsymmetry = openingAngleAsymmetry3D;
            //std::cout << "electronTreeVariables.m_postShowerStartOpeningAngleAsymmetry: " << electronTreeVariables.m_postShowerStartOpeningAngleAsymmetry << std::endl;
        }
        catch (...)
        {
        }

        const CartesianVector directionAxis(connectionPathwayDirection3D);
        CartesianVector orthoAxis1(directionAxis.GetCrossProduct(CartesianVector(1.f, 1.f, 1.f)).GetUnitVector());
        const CartesianVector orthoAxis2(directionAxis.GetCrossProduct(orthoAxis1).GetUnitVector());

        // Get angular asymmetry of hits (i don't think that there is any reason for this to be energy weighted)

        float meanTransverseAngle(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float angleToOrthoAxis1(orthoAxis1.GetOpeningAngle(projectedPosition));
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                angleToOrthoAxis1 = M_PI + (M_PI - angleToOrthoAxis1);

            meanTransverseAngle += angleToOrthoAxis1;
        }

        meanTransverseAngle /= postShowerHitList3D.size();
        meanTransverseAngle = meanTransverseAngle * 180.f / M_PI;
        electronTreeVariables.m_postShowerStartMeanTransverseAngle = meanTransverseAngle;
        //std::cout << "electronTreeVariables.m_postShowerStartMeanTransverseAngle: " << electronTreeVariables.m_postShowerStartMeanTransverseAngle << std::endl;

        float meanLWeightedTransverseAngle(0.f);
        float totalLWeight(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float angleToOrthoAxis1(orthoAxis1.GetOpeningAngle(projectedPosition));
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                angleToOrthoAxis1 = M_PI + (M_PI - angleToOrthoAxis1);

            meanLWeightedTransverseAngle += (angleToOrthoAxis1 * std::fabs(1.f / directionProjectionMagnitude));
            totalLWeight += std::fabs(1.f / directionProjectionMagnitude);
        }

        meanLWeightedTransverseAngle /= totalLWeight;
        meanLWeightedTransverseAngle = meanLWeightedTransverseAngle * 180.f / M_PI;
        electronTreeVariables.m_postShowerStartMeanLWeightedTransverseAngle = meanLWeightedTransverseAngle;
        //std::cout << "electronTreeVariables.m_postShowerStartMeanLWeightedTransverseAngle: " << electronTreeVariables.m_postShowerStartMeanLWeightedTransverseAngle << std::endl;

        // Get transverse spread - will be larger for showers and not for tracks and will be affected by any long tracks before shower

        float meanRadialDistance(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            meanRadialDistance += thisT;
        }

        ///////////////////////////////
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        ///////////////////////////////       

        meanRadialDistance /= postShowerHitList3D.size();
        electronTreeVariables.m_postShowerStartMeanRadialDistance = meanRadialDistance;
        //std::cout << "electronTreeVariables.m_postShowerStartMeanRadialDistance: " << electronTreeVariables.m_postShowerStartMeanRadialDistance << std::endl;

        float radialDistanceSigma(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            radialDistanceSigma += std::pow((thisT - meanRadialDistance), 2);
        }

        radialDistanceSigma = std::sqrt(radialDistanceSigma / postShowerHitList3D.size());
        electronTreeVariables.m_postShowerStartRadialDistanceSigma = radialDistanceSigma;
        //std::cout << "electronTreeVariables.m_postShowerStartRadialDistanceSigma: " << electronTreeVariables.m_postShowerStartRadialDistanceSigma << std::endl;

        float energyWeightedMeanRadialDistance(0.f);
        float totalWeight(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            energyWeightedMeanRadialDistance += thisT * std::fabs(pCaloHit->GetElectromagneticEnergy());
            totalWeight += std::fabs(pCaloHit->GetElectromagneticEnergy());
        }

        energyWeightedMeanRadialDistance /= totalWeight;
        electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance = energyWeightedMeanRadialDistance;
        //std::cout << "electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance: " << electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance << std::endl;

        float energyWeightedRadialDistanceSigma(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            energyWeightedRadialDistanceSigma += (std::pow(thisT - energyWeightedMeanRadialDistance, 2) * std::fabs(pCaloHit->GetElectromagneticEnergy()));
        }

        //std::cout << "energyWeightedRadialDistanceSigma: " << energyWeightedRadialDistanceSigma << std::endl;
        //std::cout << "totalWeight: " << totalWeight << std::endl;

        energyWeightedRadialDistanceSigma = std::sqrt(energyWeightedRadialDistanceSigma / totalWeight);

        const float estimatedMoliereRadius(energyWeightedRadialDistanceSigma * 1.65f);
        electronTreeVariables.m_postShowerStartEstimatedMoliereRadius = estimatedMoliereRadius;
        //std::cout << "electronTreeVariables.m_postShowerStartEstimatedMoliereRadius: " << electronTreeVariables.m_postShowerStartEstimatedMoliereRadius << std::endl;

        float lWeightedMeanRadialDistance(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            lWeightedMeanRadialDistance += thisT * std::fabs(1.f / directionProjectionMagnitude);
        }

        lWeightedMeanRadialDistance /= totalLWeight;
        electronTreeVariables.m_postShowerStartLWeightedMeanRadialDistance = lWeightedMeanRadialDistance;
        //std::cout << "electronTreeVariables.m_postShowerStartLWeightedMeanRadialDistance: " << electronTreeVariables.m_postShowerStartLWeightedMeanRadialDistance << std::endl;

        float lWeightedRadialDistanceSigma(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - nuVertexPosition);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            lWeightedRadialDistanceSigma += (std::pow(thisT - lWeightedMeanRadialDistance, 2) * std::fabs(1.f / directionProjectionMagnitude));
        }

        lWeightedRadialDistanceSigma = std::sqrt(lWeightedRadialDistanceSigma / totalLWeight);

        electronTreeVariables.m_postShowerStartLWeightedRadialDistanceSigma = lWeightedRadialDistanceSigma;
        //std::cout << "electronTreeVariables.m_postShowerStartLWeightedRadialDistanceSigma: " << electronTreeVariables.m_postShowerStartLWeightedRadialDistanceSigma << std::endl;

        // now look for any gaps (do this from the middle shower start)
        const bool isShowerDownstream((middleShowerStart3D - fullSlidingFitResult3D.GetGlobalMinLayerPosition()).GetMagnitude() < 
            (middleShowerStart3D - fullSlidingFitResult3D.GetGlobalMaxLayerPosition()).GetMagnitude());
        const CartesianVector &fitShowerDirection3D(fullSlidingFitResult3D.GetGlobalMinLayerDirection() * (isShowerDownstream ? 1.f : -1.f));

        FloatVector longitudinalProjections;

        for (const CaloHit * const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector displacement(pCaloHit->GetPositionVector() - middleShowerStart3D);
            const float longitudinalProjection(fitShowerDirection3D.GetDotProduct(displacement));

            if (longitudinalProjection > 0.f)
                longitudinalProjections.push_back(longitudinalProjection);
        }

        std::sort(longitudinalProjections.begin(), longitudinalProjections.end());

        const float initialGapSize(longitudinalProjections[0]);
        electronTreeVariables.m_postShowerStartInitialGapSize = initialGapSize;
        //std::cout << "electronTreeVariables.m_postShowerStartInitialGapSize: " << electronTreeVariables.m_postShowerStartInitialGapSize << std::endl;

        float maxGapSize(-10.f);

        for (unsigned int i = 1; i < longitudinalProjections.size(); ++i)
        {
            if ((i != 1) && (longitudinalProjections[i] > 10.f))
                continue;

            const float gapSize(longitudinalProjections[i] - longitudinalProjections[i - 1]);

            if (gapSize > maxGapSize)
                maxGapSize = gapSize;
        }

        electronTreeVariables.m_postShowerStartMaxGapSize = maxGapSize;
        //std::cout << "electronTreeVariables.m_postShowerStartMaxGapSize: " << electronTreeVariables.m_postShowerStartMaxGapSize << std::endl;
    }
    catch (...)
    {
        //std::cout << "one of the fits didnt work... :(" << std::endl;
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindPostShowerStart2DVariables(pandora::Algorithm *const pAlgorithm, const CaloHitList &spineHitList, const ParticleFlowObject *const pShowerPfo,
    const HitType hitType, const bool isDownstream3D, const CartesianVector &projectedNuVertex, const CartesianVector &projectedShowerStart, 
    const CartesianVector &initialShowerDirection, const CartesianVector &connectionPathwayDirection, int &nHits, float &scatterAngle,
    float &openingAngle, float &openingAngleAsymmetry, float &nuVertexHitAsymmetry, float &nuVertexEnergyAsymmetry, float &showerStartHitAsymmetry, 
    float &showerStartEnergyAsymmetry, float &nuVertexMeanRadialDistance, float &nuVertexEnergyWeightedMeanRadialDistance, float &showerStartMeanRadialDistance, 
    float &showerStartEnergyWeightedMeanRadialDistance, float &nuVertexMoliereRadius, float &showerStartMoliereRadius, float &positiveOpeningAngle, 
    float &negativeOpeningAngle, float &showerApexL, float &showerApexT, float &fitShowerStartL, float &fitShowerStartT, float &foundHitRatio)
{
    const bool isDownstream(projectedShowerStart.GetZ() > projectedNuVertex.GetZ());

    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, caloHitList);

    CaloHitList postShowerHitList;
    CartesianPointVector postShowerPositions;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacement(hitPosition - projectedShowerStart);
        const float l(initialShowerDirection.GetDotProduct(displacement));
        const float t(initialShowerDirection.GetCrossProduct(displacement).GetMagnitude());

        // used to be t < 14.f, should be using the opening angle?? WHO KNOWS
        if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (t < 14.f))
        {
            postShowerHitList.push_back(pCaloHit);
            postShowerPositions.push_back(pCaloHit->GetPositionVector());
        }
    }

    nHits = postShowerHitList.size();

    int foundHits(spineHitList.size());
    for (const CaloHit *const pCaloHit : postShowerHitList)
    {
        if (std::find(spineHitList.begin(), spineHitList.end(), pCaloHit) == spineHitList.end())
            ++foundHits;
    }

    foundHitRatio = static_cast<float>(foundHits) / static_cast<float>(caloHitList.size());

    std::cout << "postShowerPositions.size(): " << postShowerPositions.size() << std::endl;

    try
    {
        TwoDSlidingFitResult slidingFitResult(&postShowerPositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (std::find(postShowerHitList.begin(), postShowerHitList.end(), pCaloHit) != postShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - projectedShowerStart);
            const float l(slidingFitResult.GetGlobalMinLayerDirection().GetDotProduct(displacement));
            const float t(slidingFitResult.GetGlobalMinLayerDirection().GetCrossProduct(displacement).GetMagnitude());

            // used to be t < 14.f, should be using the opening angle?? WHO KNOWS
            if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (t < 14.f))
            {
                postShowerHitList.push_back(pCaloHit);
                postShowerPositions.push_back(pCaloHit->GetPositionVector());
            }
        }

        std::cout << "postShowerPositions.sizE(): " << postShowerPositions.size() << std::endl;

        nHits = postShowerHitList.size();

        foundHits = spineHitList.size();
        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            if (std::find(spineHitList.begin(), spineHitList.end(), pCaloHit) == spineHitList.end())
                ++foundHits;
        }

        foundHitRatio = static_cast<float>(foundHits) / static_cast<float>(caloHitList.size());

        ////////////////////////////////////////
        /*
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &postShowerHitList, "postShowerHitList", BLUE);
        const CartesianVector start1(projectedShowerStart - (initialShowerDirection * 100.f));
        const CartesianVector end1(projectedShowerStart + (initialShowerDirection * 100.f));
        const CartesianVector start2(projectedShowerStart - (slidingFitResult.GetGlobalMinLayerDirection() * 100.f));
        const CartesianVector end2(projectedShowerStart + (slidingFitResult.GetGlobalMinLayerDirection() * 100.f));
        const CartesianVector start3(slidingFitResult.GetGlobalMinLayerPosition() - (slidingFitResult.GetGlobalMinLayerDirection() * 100.f));
        const CartesianVector end3(slidingFitResult.GetGlobalMinLayerPosition() + (slidingFitResult.GetGlobalMinLayerDirection() * 100.f));
        const CartesianVector start4(projectedNuVertex - (connectionPathwayDirection*200.f));
        const CartesianVector end4(projectedNuVertex + (connectionPathwayDirection*200.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &start1, &end1, "INITIAL DIRECTION", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &start2, &end2, "FIT DIRECTION", BLACK, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &start3, &end3, "FIT DIRECTION FROM FIT POSIION", VIOLET, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &start4, &end4, "BEGINNING", BLUE, 2, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedShowerStart, "projectedShowerStart", GREEN, 2);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        ///////////////////////////////////////

        const bool isShowerDownstream((projectedShowerStart - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude() < 
            (projectedShowerStart - slidingFitResult.GetGlobalMaxLayerPosition()).GetMagnitude());

        const CartesianVector streamCorrectedConnectionPathwayDirection(connectionPathwayDirection * (isDownstream ? 1.f : -1.f));
        const CartesianVector streamCorrectedShowerDirection(isShowerDownstream ? slidingFitResult.GetGlobalMinLayerDirection() : slidingFitResult.GetGlobalMaxLayerDirection() * (-1.f));

        scatterAngle = streamCorrectedConnectionPathwayDirection.GetOpeningAngle(streamCorrectedShowerDirection) * 180.f / M_PI;

        try
        {
             // Now try to characterise the shower // will change this so dw....
            CartesianVector directionAxis1 = slidingFitResult.GetGlobalMinLayerDirection();
            CartesianVector orthoAxis1 = directionAxis1.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f));

            std::cout << "directionAxis1: " << directionAxis1 << std::endl;
            std::cout << "orthoAxis1: " << orthoAxis1 << std::endl;

            std::map<int, float> positiveEdges, negativeEdges;

            std::cout << "projectedShowerStart: " << projectedShowerStart << std::endl;

            for (const CaloHit *const pCaloHit : postShowerHitList)
            {
                const CartesianVector position(pCaloHit->GetPositionVector() - projectedShowerStart);
                const float thisT(directionAxis1.GetCrossProduct(position).GetMagnitude());
                const float thisL(directionAxis1.GetDotProduct(position));
                const float orthoL(orthoAxis1.GetDotProduct(position));

                std::map<int, float> &edgeMap(orthoL > 0.f ? positiveEdges : negativeEdges);

                const int lIndex(std::floor(thisL / 2.f));

                edgeMap[lIndex] = (edgeMap.find(lIndex) == edgeMap.end() ? thisT : std::max(edgeMap[lIndex] , thisT));;
            }

            CartesianPointVector positiveEdgePositions, negativeEdgePositions;

            for (auto &entry : positiveEdges)
                positiveEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

            for (auto &entry : negativeEdges)
                negativeEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));


            const TwoDSlidingFitResult positiveEdgeFit(&positiveEdgePositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
            const TwoDSlidingFitResult negativeEdgeFit(&negativeEdgePositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));


        std::cout << "positiveEdgePositions.size() " << positiveEdgePositions.size() << std::endl;
        std::cout << "negativeEdgePositions.size() " << negativeEdgePositions.size() << std::endl;
            std::cout << "positiveEdges.size() " << positiveEdges.size() << std::endl;
            std::cout << "negativeEdges.size() " << negativeEdges.size() << std::endl;


            const CartesianVector positiveMinLayer(positiveEdgeFit.GetGlobalMinLayerPosition());
            const CartesianVector positiveMaxLayer(positiveEdgeFit.GetGlobalMaxLayerPosition());
            const CartesianVector negativeMinLayer(negativeEdgeFit.GetGlobalMinLayerPosition());
            const CartesianVector negativeMaxLayer(negativeEdgeFit.GetGlobalMaxLayerPosition());

            const CartesianVector globalPositiveMinLayer(projectedShowerStart + (directionAxis1 * positiveMinLayer.GetZ()) + (orthoAxis1 * positiveMinLayer.GetX()));
            const CartesianVector globalPositiveMaxLayer(projectedShowerStart + (directionAxis1 * positiveMaxLayer.GetZ()) + (orthoAxis1 * positiveMaxLayer.GetX()));
            const CartesianVector globalNegativeMinLayer(projectedShowerStart + (directionAxis1 * negativeMinLayer.GetZ()) - (orthoAxis1 * negativeMinLayer.GetX()));
            const CartesianVector globalNegativeMaxLayer(projectedShowerStart + (directionAxis1 * negativeMaxLayer.GetZ()) - (orthoAxis1 * negativeMaxLayer.GetX()));

            const float positiveGradient((positiveMaxLayer.GetZ() - positiveMinLayer.GetZ()) / (positiveMaxLayer.GetX() - positiveMinLayer.GetX()));
            const float negativeGradient((negativeMaxLayer.GetZ() - negativeMinLayer.GetZ()) / (negativeMaxLayer.GetX() - negativeMinLayer.GetX()));
            const float positiveIntercept(positiveMaxLayer.GetZ() - (positiveGradient * positiveMaxLayer.GetX()));
            const float negativeIntercept(negativeMaxLayer.GetZ() - (negativeGradient * negativeMaxLayer.GetX()));

            const float collisionT((-1.f) * (positiveIntercept - negativeIntercept) / (positiveGradient + negativeGradient));
            const float collisionL((positiveGradient * collisionT) + positiveIntercept);

            const CartesianVector collisionPoint(projectedShowerStart + (directionAxis1 * collisionL) + (orthoAxis1 * collisionT));

            const CartesianVector positiveEdgeVector((globalPositiveMaxLayer - globalPositiveMinLayer).GetUnitVector());
            const CartesianVector negativeEdgeVector((globalNegativeMaxLayer - globalNegativeMinLayer).GetUnitVector());

            positiveOpeningAngle = directionAxis1.GetOpeningAngle(positiveEdgeVector) * 180.f / M_PI;
            negativeOpeningAngle = directionAxis1.GetOpeningAngle(negativeEdgeVector) * 180.f / M_PI;
            openingAngle = positiveOpeningAngle + negativeOpeningAngle;

            const CartesianVector fitShowerStart(isShowerDownstream ? slidingFitResult.GetGlobalMinLayerPosition() : slidingFitResult.GetGlobalMaxLayerPosition());
            showerApexL = std::min(std::max(streamCorrectedShowerDirection.GetDotProduct(collisionPoint - fitShowerStart), -1000.f), 1000.f);
            showerApexT = streamCorrectedShowerDirection.GetCrossProduct(collisionPoint - fitShowerStart).GetMagnitude();

            //////////////////////////////////////
            /*
            std::cout << "showerApexL: " << showerApexL << std::endl;
            std::cout << "showerApexT: " << showerApexT << std::endl;
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &collisionPoint, "collisionPoint", BLACK, 2);
            PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &globalPositiveMinLayer, &globalPositiveMaxLayer, "POSITIVE SHOWER EDGE", BLUE, 2, 2);
            PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &globalNegativeMinLayer, &globalNegativeMaxLayer, "NEGATIVE SHOWER EDGE", BLUE, 2, 2);
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            */
            //////////////////////////////////////
        }
        catch(...)
        {
            //std::cout << "CAUGHT" << std::endl;
            positiveOpeningAngle = -10.f;
            negativeOpeningAngle = -10.f;
            openingAngle = -10.f;
            showerApexL = -9999.f;
            showerApexT = -10.f;
            foundHitRatio = -0.5f;
        }

        CartesianVector directionAxis1 = connectionPathwayDirection;
        CartesianVector orthoAxis1 = directionAxis1.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f));

        // does it look like the shower is coming from the nu vertex?
        const CartesianVector fitShowerStart(isShowerDownstream ? slidingFitResult.GetGlobalMinLayerPosition() : slidingFitResult.GetGlobalMaxLayerPosition());

        fitShowerStartL = directionAxis1.GetDotProduct(fitShowerStart - projectedNuVertex);
        fitShowerStartT = directionAxis1.GetCrossProduct(fitShowerStart - projectedNuVertex).GetMagnitude();

        CartesianVector directionAxis2 = streamCorrectedShowerDirection;
        CartesianVector orthoAxis2 = directionAxis2.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f));

        nuVertexHitAsymmetry = 0.f;

        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - projectedNuVertex);
            const float thisL(orthoAxis1.GetDotProduct(position));

            nuVertexHitAsymmetry += (thisL < 0.f) ? -1.f : 1.f;
        }

        nuVertexHitAsymmetry = std::fabs(nuVertexHitAsymmetry);

        float totalEnergy(0.f);
        nuVertexEnergyAsymmetry = 0.f;

        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));

            totalEnergy += hitEnergy;

            const CartesianVector position(pCaloHit->GetPositionVector() - projectedNuVertex);
            const float thisL(orthoAxis1.GetDotProduct(position));

            nuVertexEnergyAsymmetry += (thisL < 0.f) ? (-1.f * hitEnergy) : hitEnergy;
        }

        nuVertexEnergyAsymmetry = (totalEnergy < std::numeric_limits<float>::epsilon()) ? -0.5f : (nuVertexEnergyAsymmetry / totalEnergy);
        nuVertexEnergyAsymmetry = std::fabs(nuVertexEnergyAsymmetry);

        showerStartHitAsymmetry = 0.f;

        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fitShowerStart);
            const float thisL(orthoAxis2.GetDotProduct(position));

            showerStartHitAsymmetry += (thisL < 0.f) ? -1.f : 1.f;
        }

        showerStartHitAsymmetry = std::fabs(showerStartHitAsymmetry);

        showerStartEnergyAsymmetry = 0.f;

        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));
            const CartesianVector position(pCaloHit->GetPositionVector() - fitShowerStart);
            const float thisL(orthoAxis2.GetDotProduct(position));

            showerStartEnergyAsymmetry += (thisL < 0.f) ? (-1.f * hitEnergy) : hitEnergy;
        }

        showerStartEnergyAsymmetry = (totalEnergy < std::numeric_limits<float>::epsilon()) ? -0.5f : (showerStartEnergyAsymmetry / totalEnergy);
        showerStartEnergyAsymmetry = std::fabs(showerStartEnergyAsymmetry);

        // Get mean radial distance

        nuVertexMeanRadialDistance = 0.f;
        nuVertexEnergyWeightedMeanRadialDistance = 0.f;
        showerStartMeanRadialDistance = 0.f;
        showerStartEnergyWeightedMeanRadialDistance = 0.f;

        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            const CartesianVector nuPosition(pCaloHit->GetPositionVector() - projectedNuVertex);
            const CartesianVector showerStartPosition(pCaloHit->GetPositionVector() - fitShowerStart);
            const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));

            nuVertexMeanRadialDistance += directionAxis1.GetCrossProduct(nuPosition).GetMagnitude();
            nuVertexEnergyWeightedMeanRadialDistance += (directionAxis1.GetCrossProduct(nuPosition).GetMagnitude() * hitEnergy);
            showerStartMeanRadialDistance += directionAxis2.GetCrossProduct(showerStartPosition).GetMagnitude();
            showerStartEnergyWeightedMeanRadialDistance += (directionAxis2.GetCrossProduct(showerStartPosition).GetMagnitude() * hitEnergy);
        }

        nuVertexMeanRadialDistance = (postShowerHitList.size() > 0) ? nuVertexMeanRadialDistance / postShowerHitList.size() : -10.f;
        nuVertexEnergyWeightedMeanRadialDistance = (totalEnergy < std::numeric_limits<float>::epsilon()) ? -10.f : nuVertexEnergyWeightedMeanRadialDistance / totalEnergy;
        showerStartMeanRadialDistance = (postShowerHitList.size() > 0) ? showerStartMeanRadialDistance / postShowerHitList.size() : -10.f;
        showerStartEnergyWeightedMeanRadialDistance = (totalEnergy < std::numeric_limits<float>::epsilon()) ? -10.f : showerStartEnergyWeightedMeanRadialDistance / totalEnergy;

        //std::cout << "111111111111" << std::endl;
        // Get Molliere radius...

        CaloHitVector nuVertexPostShowerHitVector(postShowerHitList.begin(), postShowerHitList.end());

        std::sort(nuVertexPostShowerHitVector.begin(), nuVertexPostShowerHitVector.end(),
            [&projectedNuVertex, &directionAxis1](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool {
                const CartesianVector positionA(pCaloHitA->GetPositionVector() - projectedNuVertex);
                const CartesianVector positionB(pCaloHitB->GetPositionVector() - projectedNuVertex);

                const float tA(directionAxis1.GetCrossProduct(positionA).GetMagnitude());
                const float tB(directionAxis1.GetCrossProduct(positionB).GetMagnitude());

                return tA < tB;
            });

        float nuVertexRunningEnergySum(0.f);
        nuVertexMoliereRadius = -10.f;

        for (const CaloHit *const pCaloHit : nuVertexPostShowerHitVector)
        {
            const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));
            nuVertexRunningEnergySum += hitEnergy;

            if ((nuVertexRunningEnergySum / totalEnergy) > 0.9f)
            {
                const CartesianVector position(pCaloHit->GetPositionVector() - projectedNuVertex);
                nuVertexMoliereRadius = directionAxis1.GetCrossProduct(position).GetMagnitude();
                break;
            }
        }

        CaloHitVector showerStartPostShowerHitVector(postShowerHitList.begin(), postShowerHitList.end());

        std::sort(showerStartPostShowerHitVector.begin(), showerStartPostShowerHitVector.end(),
            [&fitShowerStart, &directionAxis2](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool {
                const CartesianVector positionA(pCaloHitA->GetPositionVector() - fitShowerStart);
                const CartesianVector positionB(pCaloHitB->GetPositionVector() - fitShowerStart);

                const float tA(directionAxis2.GetCrossProduct(positionA).GetMagnitude());
                const float tB(directionAxis2.GetCrossProduct(positionB).GetMagnitude());

                return tA < tB;
            });

        float showerStartRunningEnergySum(0.f);
        showerStartMoliereRadius = -10.f;

        for (const CaloHit *const pCaloHit : showerStartPostShowerHitVector)
        {
            const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));
            showerStartRunningEnergySum += hitEnergy;

            if ((showerStartRunningEnergySum / totalEnergy) > 0.9f)
            {
                const CartesianVector position(pCaloHit->GetPositionVector() - fitShowerStart);
                showerStartMoliereRadius = directionAxis2.GetCrossProduct(position).GetMagnitude();
                break;
            }
        }
    }
    catch (...)
    {
        scatterAngle = -10.f;
        openingAngle = -10.f;
        openingAngleAsymmetry = -10.f;
        nuVertexHitAsymmetry = -10.f;
        nuVertexEnergyAsymmetry = -0.5f;
        showerStartHitAsymmetry = -10.f;
        showerStartEnergyAsymmetry = -0.5f;
        nuVertexMeanRadialDistance = -10.f;
        nuVertexEnergyWeightedMeanRadialDistance = -10.f;
        showerStartMeanRadialDistance = -10.f;
        showerStartEnergyWeightedMeanRadialDistance = -10.f;

        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillConnectionPathwayVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition, CartesianVector &minShowerStart3D,
    CartesianVector &middleShowerStart3D, CartesianVector &maxShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    std::cout << "minShowerStart3D: " << minShowerStart3D << std::endl;
    std::cout << "middleShowerStart3D: " << middleShowerStart3D << std::endl;
    std::cout << "maxShowerStart3D: " << maxShowerStart3D << std::endl;


    electronTreeVariables.m_pathwayLengthMin = std::min((nuVertexPosition - minShowerStart3D).GetMagnitude(), 30.f);
    electronTreeVariables.m_pathwayLengthMiddle = (nuVertexPosition - middleShowerStart3D).GetMagnitude();
    electronTreeVariables.m_pathwayLengthMax = (nuVertexPosition - maxShowerStart3D).GetMagnitude();
    electronTreeVariables.m_pathwayShowerStartDelta = (maxShowerStart3D - minShowerStart3D).GetMagnitude();

    try
    {
        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
            spinePositionsU.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
            spinePositionsV.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
            spinePositionsW.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        electronTreeVariables.m_pathwayMaxScatteringAngle = std::min(LArConnectionPathwayHelper::GetLargest3DKink(pAlgorithm, spineFitU, spineFitV, spineFitW,
            nuVertexPosition, maxShowerStart3D), 20.f);

        electronTreeVariables.m_minShowerStartPathwayScatteringAngle2D = LArConnectionPathwayHelper::Get2DKink(pAlgorithm, spineFitU, spineFitV, spineFitW, minShowerStart3D);
        electronTreeVariables.m_middleShowerStartPathwayScatteringAngle2D = LArConnectionPathwayHelper::Get2DKink(pAlgorithm, spineFitU, spineFitV, spineFitW, middleShowerStart3D);
        electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D = std::min(LArConnectionPathwayHelper::Get2DKink(pAlgorithm, spineFitU, spineFitV, spineFitW, maxShowerStart3D), 10.f);
    }
    catch(...)
    {
    }

    // middle
    LArConnectionPathwayHelper::FillConnectionPathwayEnergyVariables(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, minShowerStart3D, electronTreeVariables);


    // here, think about moving the shower vertex... 
    /*
    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    try
    {
        // try 10 so i can make it consistent?
        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const CartesianVector &startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection());
        const CartesianVector &startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection());
        const CartesianVector &startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection());

        const CartesianVector projectedShowerStartU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_U));
        const CartesianVector projectedShowerStartV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_V));
        const CartesianVector projectedShowerStartW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_W));

        CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
            protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
        {
            if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
                protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, protoShowerU.m_connectionPathway.m_startDirection, protoShowerV.m_connectionPathway.m_startDirection, 
                protoShowerW.m_connectionPathway.m_startDirection, connectionPathwayDirection3D))
            {
                //std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
                return false;
            }
        }

        CaloHitList caloHitList3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHitList3D);

        CaloHitList postShowerHitList3D; 
        CartesianPointVector postShowerPositions3D;

        float highestL(0.f);
        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - minShowerStart3D);
            const float l(connectionPathwayDirection3D.GetDotProduct(displacement));
            const float t(connectionPathwayDirection3D.GetCrossProduct(displacement).GetMagnitude());

            const float openingAngle(displacement.GetOpeningAngle(connectionPathwayDirection3D) * 180.f / M_PI);

            if ((l > 0.f) && (l > highestL) && (openingAngle < 25.f) && (t < 10.f))
                highestL = l;
        }

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - minShowerStart3D);
            const float l(connectionPathwayDirection3D.GetDotProduct(displacement));
            const float t(connectionPathwayDirection3D.GetCrossProduct(displacement).GetMagnitude());

            const float openingAngle(displacement.GetOpeningAngle(connectionPathwayDirection3D) * 180.f / M_PI);

            if ((l > 0.f) && (l < (0.4 * highestL)) && (openingAngle < 25.f) && (t < 10.f))
            {
                postShowerHitList3D.push_back(pCaloHit);
                postShowerPositions3D.push_back(pCaloHit->GetPositionVector());
            }
        }

        const CartesianVector start4(minShowerStart3D - (connectionPathwayDirection3D * 1000.f));
        const CartesianVector end4(minShowerStart3D + (connectionPathwayDirection3D * 1000.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &start4, &end4, "DIRECTION AXIS", VIOLET, 2, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &minShowerStart3D, "minShowerStart3D", BLACK, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &postShowerHitList3D, "postShowerHitList3D", RED);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        ThreeDSlidingFitResult slidingFitResult3D(&postShowerPositions3D, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        CartesianVector projectedShowerStart(slidingFitResult3D.GetGlobalMinLayerPosition());
        const CartesianVector &jam(slidingFitResult3D.GetGlobalMinLayerDirection());


        highestL = 0.f;
        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - projectedShowerStart);
            const float l(jam.GetDotProduct(displacement));
            const float t(jam.GetDotProduct(displacement));

            const float openingAngle(displacement.GetOpeningAngle(jam) * 180.f / M_PI);

            if ((l > 0.f) && (l > highestL) && (openingAngle < 35.f) && (t < 10.f))
                highestL = l;
        }

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            if (std::find(postShowerHitList3D.begin(), postShowerHitList3D.end(), pCaloHit) != postShowerHitList3D.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - projectedShowerStart);
            const float l(jam.GetDotProduct(displacement));
            const float t(jam.GetDotProduct(displacement));

            const float openingAngle(displacement.GetOpeningAngle(jam) * 180.f / M_PI);

            if ((l > 0.f) && (l < (0.4*highestL)) && (openingAngle < 35.f) && (t < 10.f))
            {
                postShowerHitList3D.push_back(pCaloHit);
                postShowerPositions3D.push_back(pCaloHit->GetPositionVector());
            }
        }

        const CartesianVector start5(projectedShowerStart - (jam * 1000.f));
        const CartesianVector end5(projectedShowerStart + (jam * 1000.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &start5, &end5, "DIRECTION AXIS", VIOLET, 2, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedShowerStart, "projectedShowerStart", BLACK, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &postShowerHitList3D, "postShowerHitList3D", RED);


        ThreeDSlidingFitResult fullSlidingFitResult3D(&postShowerPositions3D, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        projectedShowerStart = fullSlidingFitResult3D.GetGlobalMinLayerPosition();

        const CartesianVector directionAxis(fullSlidingFitResult3D.GetGlobalMinLayerDirection());
        CartesianVector orthoAxis1(directionAxis.GetCrossProduct(CartesianVector(1.f, 1.f, 1.f)).GetUnitVector());
        const CartesianVector orthoAxis2(directionAxis.GetCrossProduct(orthoAxis1).GetUnitVector());

        std::map<int, float> positiveEdges, negativeEdges;

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - projectedShowerStart);
            const float thisT(directionAxis.GetCrossProduct(position).GetMagnitude());
            const float thisL(directionAxis.GetDotProduct(position));
            const float orthoL(orthoAxis1.GetDotProduct(position));

            std::map<int, float> &edgeMap(orthoL > 0.f ? positiveEdges : negativeEdges);

            const int lIndex(std::floor(thisL / 2.f));
            edgeMap[lIndex] = (edgeMap.find(lIndex) == edgeMap.end() ? thisT : std::max(edgeMap[lIndex], thisT));
        }

        CartesianPointVector positiveEdgePositions, negativeEdgePositions;

        for (auto &entry : positiveEdges) 
            positiveEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

        for (auto &entry : negativeEdges) 
            negativeEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

        const TwoDSlidingFitResult positiveEdgeFit(&positiveEdgePositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult negativeEdgeFit(&negativeEdgePositions, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const CartesianVector positiveMinLayer(positiveEdgeFit.GetGlobalMinLayerPosition());
        const CartesianVector positiveMaxLayer(positiveEdgeFit.GetGlobalMaxLayerPosition());
        const CartesianVector negativeMinLayer(negativeEdgeFit.GetGlobalMinLayerPosition());
        const CartesianVector negativeMaxLayer(negativeEdgeFit.GetGlobalMaxLayerPosition());

        const float positiveGradient((positiveMaxLayer.GetZ() - positiveMinLayer.GetZ()) / (positiveMaxLayer.GetX() - positiveMinLayer.GetX()));
        const float negativeGradient((negativeMaxLayer.GetZ() - negativeMinLayer.GetZ()) / (negativeMaxLayer.GetX() - negativeMinLayer.GetX()));
        const float positiveIntercept(positiveMaxLayer.GetZ() - (positiveGradient * positiveMaxLayer.GetX()));
        const float negativeIntercept(negativeMaxLayer.GetZ() - (negativeGradient * negativeMaxLayer.GetX()));

        std::cout << "positiveGradient: " << positiveGradient << std::endl;
        std::cout << "negativeGradient: " << negativeGradient << std::endl;
        std::cout << "positiveIntercept: " << positiveIntercept << std::endl;
        std::cout << "negativeIntercept: " << negativeIntercept << std::endl;

        const float collisionT((-1.f) * (positiveIntercept + negativeIntercept) / (positiveGradient + negativeGradient));
        const float collisionL((positiveGradient * collisionT) + positiveIntercept);

        const float positiveStartL(positiveIntercept);
        const float negativeStartL(negativeIntercept);

        std::cout << "positiveStartL: " << positiveStartL << std::endl;
        std::cout << "negativeStartL: " << negativeStartL << std::endl;

        const CartesianVector globalPositiveMinLayer(projectedShowerStart + (directionAxis * positiveMinLayer.GetZ()) + (orthoAxis1 * positiveMinLayer.GetX()));
        const CartesianVector globalPositiveMaxLayer(projectedShowerStart + (directionAxis * positiveMaxLayer.GetZ()) + (orthoAxis1 * positiveMaxLayer.GetX()));
        const CartesianVector globalNegativeMinLayer(projectedShowerStart + (directionAxis * negativeMinLayer.GetZ()) + (orthoAxis1 * negativeMinLayer.GetX()));
        const CartesianVector globalNegativeMaxLayer(projectedShowerStart + (directionAxis * negativeMaxLayer.GetZ()) + (orthoAxis1 * negativeMaxLayer.GetX()));
        const CartesianVector collisionPoint(projectedShowerStart + (directionAxis * collisionL));
        const CartesianVector positiveStart(projectedShowerStart + (directionAxis * positiveStartL));
        const CartesianVector negativeStart(projectedShowerStart + (directionAxis * negativeStartL));

        const CartesianVector positiveEdgeVector((globalPositiveMaxLayer - globalPositiveMinLayer).GetUnitVector());
        const CartesianVector negativeEdgeVector((globalNegativeMaxLayer - globalNegativeMinLayer).GetUnitVector());

        const float positiveOpeningAngle(directionAxis.GetOpeningAngle(positiveEdgeVector));
        const float negativeOpeningAngle(directionAxis.GetOpeningAngle(negativeEdgeVector));

        std::cout << "positiveOpeningAngle: " << positiveOpeningAngle << std::endl;
        std::cout << "negativeOpeningAngle: " << negativeOpeningAngle << std::endl;

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveStart, "positiveStart", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &collisionPoint, "collisionPoint", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeStart, "negativeStart", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        
        // Get angular asymmetry of hits (i don't think that there is any reason for this to be energy weighted)

        float meanTransverseAngle(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - minShowerStart3D);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float angleToOrthoAxis1(orthoAxis1.GetOpeningAngle(projectedPosition));
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                angleToOrthoAxis1 = M_PI + (M_PI - angleToOrthoAxis1);

            meanTransverseAngle += angleToOrthoAxis1;
        }

        meanTransverseAngle /= postShowerHitList3D.size();
        meanTransverseAngle = meanTransverseAngle * 180.f / M_PI;
        std::cout << "meanTransverseAngle: " << meanTransverseAngle << std::endl;

        float energyWeightedMeanRadialDistance(0.f);
        float totalWeight(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - minShowerStart3D);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            energyWeightedMeanRadialDistance += thisT * std::fabs(pCaloHit->GetElectromagneticEnergy());
            totalWeight += std::fabs(pCaloHit->GetElectromagneticEnergy());
        }

        energyWeightedMeanRadialDistance /= totalWeight;
        std::cout << "energyWeightedMeanRadialDistance: " << energyWeightedMeanRadialDistance << std::endl;

        float energyWeightedRadialDistanceSigma(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - minShowerStart3D);
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            energyWeightedRadialDistanceSigma += (std::pow(thisT - energyWeightedMeanRadialDistance, 2) * std::fabs(pCaloHit->GetElectromagneticEnergy()));
        }

        energyWeightedRadialDistanceSigma = std::sqrt(energyWeightedRadialDistanceSigma / totalWeight);

        const float estimatedMoliereRadius(energyWeightedRadialDistanceSigma * 1.65f);
        std::cout << "estimatedMoliereRadius: " << estimatedMoliereRadius << std::endl;
    }
    catch(...)
    {
        std::cout << "izzie it didnt work :( " << std::endl;
    }
    */

    //////////////////////////////
    /*    
    std::cout << "electronTreeVariables.m_pathwayLengthMin: " << electronTreeVariables.m_pathwayLengthMin << std::endl;
    std::cout << "electronTreeVariables.m_pathwayLengthMiddle: " << electronTreeVariables.m_pathwayLengthMiddle << std::endl;
    std::cout << "electronTreeVariables.m_pathwayLengthMax: " << electronTreeVariables.m_pathwayLengthMax << std::endl;
    std::cout << "electronTreeVariables.m_pathwayShowerStartDelta: " << electronTreeVariables.m_pathwayShowerStartDelta << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMaxScatteringAngle: " << electronTreeVariables.m_pathwayMaxScatteringAngle << std::endl;
    std::cout << "electronTreeVariables.m_minShowerStartPathwayScatteringAngle2D: " << electronTreeVariables.m_minShowerStartPathwayScatteringAngle2D << std::endl;
    std::cout << "electronTreeVariables.m_middleShowerStartPathwayScatteringAngle2D: " << electronTreeVariables.m_middleShowerStartPathwayScatteringAngle2D << std::endl;
    std::cout << "electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D: " << electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergyMeanU: " << electronTreeVariables.m_pathwayEnergyMeanU << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergyMeanV: " << electronTreeVariables.m_pathwayEnergyMeanV << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergyMeanW: " << electronTreeVariables.m_pathwayEnergyMeanW << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergySigmaU: " << electronTreeVariables.m_pathwayEnergySigmaU << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergySigmaV: " << electronTreeVariables.m_pathwayEnergySigmaV << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergySigmaW: " << electronTreeVariables.m_pathwayEnergySigmaW << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMinEnergyMean: " << electronTreeVariables.m_pathwayMinEnergyMean << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMiddleEnergyMean: " << electronTreeVariables.m_pathwayMiddleEnergyMean << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMaxEnergyMean: " << electronTreeVariables.m_pathwayMaxEnergyMean << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMinEnergySigma: " << electronTreeVariables.m_pathwayMinEnergySigma << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMiddleEnergySigma: " << electronTreeVariables.m_pathwayMiddleEnergySigma << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMaxEnergySigma: " << electronTreeVariables.m_pathwayMaxEnergySigma << std::endl;
    */
    /////////////////////////////

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::GetLargest3DKink(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFitU, const TwoDSlidingFitResult &spineFitV, 
    const TwoDSlidingFitResult &spineFitW, const CartesianVector &nuVertexPosition, CartesianVector &maxShowerStart3D)
{
    float maxOpeningAngle(-std::numeric_limits<float>::max());

    maxOpeningAngle = std::max(maxOpeningAngle, LArConnectionPathwayHelper::GetLargest3DKinkFromView(pAlgorithm, spineFitU, spineFitV, spineFitW, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, maxShowerStart3D));
    maxOpeningAngle = std::max(maxOpeningAngle, LArConnectionPathwayHelper::GetLargest3DKinkFromView(pAlgorithm, spineFitV, spineFitW, spineFitU, TPC_VIEW_V, TPC_VIEW_W, TPC_VIEW_U, maxShowerStart3D));
    maxOpeningAngle = std::max(maxOpeningAngle, LArConnectionPathwayHelper::GetLargest3DKinkFromView(pAlgorithm, spineFitW, spineFitU, spineFitV, TPC_VIEW_W, TPC_VIEW_U, TPC_VIEW_V, maxShowerStart3D));

    return maxOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------


float LArConnectionPathwayHelper::GetLargest3DKinkFromView(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, 
    const TwoDSlidingFitResult &spineFit1, const TwoDSlidingFitResult &spineFit2, const HitType hitType, const HitType hitType1, 
    const HitType hitType2, const CartesianVector &maxShowerStart3D)
{
    const LayerFitResultMap &layerFitResultMap(spineFit.GetLayerFitResultMap());
    const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);

    const int nLayersHalfWindow(spineFit.GetLayerFitHalfWindow());
    const int nLayersSpanned(1 + maxLayer - minLayer);

    //std::cout << "minLayer: " << minLayer << std::endl;
    //std::cout << "maxLayer: " << maxLayer << std::endl;
    //std::cout << "nLayersHalfWindow: " << nLayersHalfWindow << std::endl;

    if (nLayersSpanned <= 2 * nLayersHalfWindow)
    {
        //std::cout << "failed this window cut" << std::endl;
        return -10.f;
    }

    const CartesianVector projectedShowerStart(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), maxShowerStart3D, hitType));

    float showerStartL(0.f), showerStartT(0.f);
    spineFit.GetLocalPosition(projectedShowerStart, showerStartL, showerStartT);
    float maxCentralLayer(spineFit.GetLayer(showerStartL) - nLayersHalfWindow);

    float highestOpeningAngle(-10.f);
    CartesianVector kinkPosition(0.f, 0.f, 0.f);

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        const int iLayer(iter->first);

        if (iLayer > maxCentralLayer)
            continue;

        const float rL(spineFit.GetL(iLayer));
        const float rL1(spineFit.GetL(iLayer - nLayersHalfWindow));
        const float rL2(spineFit.GetL(iLayer + nLayersHalfWindow));

        CartesianVector firstPosition(0.f, 0.f, 0.f), centralPosition(0.f, 0.f, 0.f), secondPosition(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL1, firstPosition)) ||
            (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL, centralPosition)) ||
            (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL2, secondPosition)))
       {
           continue;
       }

        float firstPositionX(firstPosition.GetX()), centralPositionX(centralPosition.GetX()), secondPositionX(secondPosition.GetX());

        CartesianVector firstPosition1(0.f, 0.f, 0.f), centralPosition1(0.f, 0.f, 0.f), secondPosition1(0.f, 0.f, 0.f);
        CartesianVector firstPosition2(0.f, 0.f, 0.f), centralPosition2(0.f, 0.f, 0.f), secondPosition2(0.f, 0.f, 0.f);

        if ((spineFit1.GetGlobalFitPositionAtX(firstPositionX, firstPosition1) != STATUS_CODE_SUCCESS) ||
            (spineFit1.GetGlobalFitPositionAtX(centralPositionX, centralPosition1) != STATUS_CODE_SUCCESS) ||
            (spineFit1.GetGlobalFitPositionAtX(secondPositionX, secondPosition1) != STATUS_CODE_SUCCESS) ||
            (spineFit2.GetGlobalFitPositionAtX(firstPositionX, firstPosition2) != STATUS_CODE_SUCCESS) ||
            (spineFit2.GetGlobalFitPositionAtX(centralPositionX, centralPosition2) != STATUS_CODE_SUCCESS) ||
            (spineFit2.GetGlobalFitPositionAtX(secondPositionX, secondPosition2) != STATUS_CODE_SUCCESS))
        {
            continue;
        }

        float metric(0.f);

        const CartesianVector &firstPositionU(hitType == TPC_VIEW_U ? firstPosition : hitType1 == TPC_VIEW_U ? firstPosition1 : firstPosition2);
        const CartesianVector &centralPositionU(hitType == TPC_VIEW_U ? centralPosition : hitType1 == TPC_VIEW_U ? centralPosition1 : centralPosition2);
        const CartesianVector &secondPositionU(hitType == TPC_VIEW_U ? secondPosition : hitType1 == TPC_VIEW_U ? secondPosition1 : secondPosition2);

        /////////////////////////////////////////
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &firstPositionU, "firstPositionU", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &centralPositionU, "centralPositionU", ORANGE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &secondPositionU, "secondPositionU", GREEN, 2);
        */
        /////////////////////////////////////////

        const CartesianVector &firstPositionV(hitType == TPC_VIEW_V ? firstPosition : hitType1 == TPC_VIEW_V ? firstPosition1 : firstPosition2);
        const CartesianVector &centralPositionV(hitType == TPC_VIEW_V ? centralPosition : hitType1 == TPC_VIEW_V ? centralPosition1 : centralPosition2);
        const CartesianVector &secondPositionV(hitType == TPC_VIEW_V ? secondPosition : hitType1 == TPC_VIEW_V ? secondPosition1 : secondPosition2);

        /////////////////////////////////////////
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &firstPositionV, "firstPositionV", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &centralPositionV, "centralPositionV", ORANGE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &secondPositionV, "secondPositionV", GREEN, 2);
        */
        /////////////////////////////////////////

        const CartesianVector &firstPositionW(hitType == TPC_VIEW_W ? firstPosition : hitType1 == TPC_VIEW_W ? firstPosition1 : firstPosition2);
        const CartesianVector &centralPositionW(hitType == TPC_VIEW_W ? centralPosition : hitType1 == TPC_VIEW_W ? centralPosition1 : centralPosition2);
        const CartesianVector &secondPositionW(hitType == TPC_VIEW_W ? secondPosition : hitType1 == TPC_VIEW_W ? secondPosition1 : secondPosition2);

        /////////////////////////////////////////
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &firstPositionW, "firstPositionW", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &centralPositionW, "centralPositionW", ORANGE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &secondPositionW, "secondPositionW", GREEN, 2);
        */
        /////////////////////////////////////////

        if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, firstPositionU, firstPositionV, firstPositionW, 5.f, 2.f, metric))
        {
            //std::cout << "NOT CONSISTENT" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            continue;
        }

        if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, centralPositionU, centralPositionV, centralPositionW, 5.f, 2.f, metric))
        {
            //std::cout << "NOT CONSISTENT" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            continue;
        }

        if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, secondPositionU, secondPositionV, secondPositionW, 5.f, 2.f, metric))
        {
            //std::cout << "NOT CONSISTENT" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            continue;
        }

        float chi2(0.0);
        CartesianVector firstPosition3D(0.f, 0.f, 0.f), centralPosition3D(0.f, 0.f, 0.f), secondPosition3D(0.f, 0.f, 0.f);

        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, firstPositionU, firstPositionV, firstPositionW, firstPosition3D, chi2);
        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, centralPositionU, centralPositionV, centralPositionW, centralPosition3D, chi2);
        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, secondPositionU, secondPositionV, secondPositionW, secondPosition3D, chi2);

        const CartesianVector firstDirection3D((centralPosition3D - firstPosition3D).GetUnitVector());
        const CartesianVector secondDirection3D((secondPosition3D - centralPosition3D).GetUnitVector());

        const float openingAngle3D(secondDirection3D.GetOpeningAngle(firstDirection3D) * 180.f / M_PI);

        //std::cout << "openingAngle3D: " << openingAngle3D << std::endl;
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        if (openingAngle3D > highestOpeningAngle)
        {
            highestOpeningAngle = openingAngle3D;
            kinkPosition = centralPosition3D;
        }
    }

    /////////////////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &kinkPosition, "kinkPosition", RED, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////////////

    return highestOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::Get2DKink(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFitU, 
    const TwoDSlidingFitResult &spineFitV, const TwoDSlidingFitResult &spineFitW, const CartesianVector &showerStart3D)
{
    const float scatterAngleU(LArConnectionPathwayHelper::GetLargest2DKinkFromView(pAlgorithm, spineFitU, TPC_VIEW_U, showerStart3D));
    const float scatterAngleV(LArConnectionPathwayHelper::GetLargest2DKinkFromView(pAlgorithm, spineFitV, TPC_VIEW_V, showerStart3D));
    const float scatterAngleW(LArConnectionPathwayHelper::GetLargest2DKinkFromView(pAlgorithm, spineFitW, TPC_VIEW_W, showerStart3D));

    const float minScatterAngle(std::min(std::min(scatterAngleU, scatterAngleV), scatterAngleW));
    const float maxScatterAngle(std::max(std::max(scatterAngleU, scatterAngleV), scatterAngleW));

    for (const float scatterAngle : {scatterAngleU, scatterAngleV, scatterAngleW})
    {
        if ((std::fabs(scatterAngle - maxScatterAngle) < std::numeric_limits<float>::epsilon()) || 
            (std::fabs(scatterAngle - minScatterAngle) < std::numeric_limits<float>::epsilon()))
        {
            continue;
        }

        return scatterAngle;
    }

    return -10.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::GetLargest2DKinkFromView(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, 
    const HitType hitType, const CartesianVector &showerStart3D)
{
    const LayerFitResultMap &layerFitResultMap(spineFit.GetLayerFitResultMap());
    const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);

    const int nLayersHalfWindow(5);
    const int nLayersSpanned(1 + maxLayer - minLayer);

    if (nLayersSpanned <= 2 * nLayersHalfWindow)
    {
        std::cout << "NOT ENOUGH LAYERS SPANNED" << std::endl;
        return -10.f;
    }

    const CartesianVector projectedShowerStart(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), showerStart3D, hitType));

    float showerStartL(0.f), showerStartT(0.f);
    spineFit.GetLocalPosition(projectedShowerStart, showerStartL, showerStartT);

    float maxCentralLayer(spineFit.GetLayer(showerStartL) - nLayersHalfWindow);
    float minCentralLayer(layerFitResultMap.begin()->first + nLayersHalfWindow + 1);

    float highestOpeningAngle(-10.f);
    CartesianVector kinkPosition(0.f, 0.f, 0.f);
    CartesianVector highestKinkPosition(0.f, 0.f, 0.f);

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        const int layer(iter->first);

        if (layer < minCentralLayer)
            continue;

        if (layer > maxCentralLayer)
            continue;

        bool found(false);
        float openingAngle2D(std::numeric_limits<float>::max());

        float thisHighestOpeningAngle(-10.f);
        CartesianVector thisHighestKinkPosition(0.f, 0.f, 0.f);
        for (int i = 0; i < nLayersHalfWindow; ++i)
        {
            const int testLayer(layer + i);
            const float rL(spineFit.GetL(testLayer));
            const float rL1(spineFit.GetL(testLayer - nLayersHalfWindow));
            const float rL2(spineFit.GetL(testLayer + nLayersHalfWindow));

            CartesianVector firstPosition(0.f, 0.f, 0.f), centralPosition(0.f, 0.f, 0.f), secondPosition(0.f, 0.f, 0.f);

            if ((STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL1, firstPosition)) ||
                (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL, centralPosition)) ||
                (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL2, secondPosition)))
            {
                continue;
            }

            const CartesianVector firstDirection((centralPosition - firstPosition).GetUnitVector());
            const CartesianVector secondDirection((secondPosition - centralPosition).GetUnitVector());

            const float testOpeningAngle2D(secondDirection.GetOpeningAngle(firstDirection) * 180.f / M_PI);

            //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &centralPosition, "centralPosition", ORANGE, 2);

            if (testOpeningAngle2D < openingAngle2D)
            {
                openingAngle2D = testOpeningAngle2D;
                found = true;
            }

            if (testOpeningAngle2D > thisHighestOpeningAngle)
            {
                thisHighestOpeningAngle = openingAngle2D;
                thisHighestKinkPosition = centralPosition;
            }
        }

        if (found)
        {
            ////////////////////////////
            //std::cout << "openingAngle2D: " << openingAngle2D << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ////////////////////////////
        }

        if (found)
        {
            if (openingAngle2D > highestOpeningAngle)
            {
                highestOpeningAngle = std::max(highestOpeningAngle, openingAngle2D);
                highestKinkPosition = thisHighestKinkPosition;
            }
        }
    }

    ///////////////////////////////
    /*
    std::cout << "highestOpeningAngle: " << highestOpeningAngle << std::endl;
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &highestKinkPosition, "highestKinkPosition", RED, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    ///////////////////////////////

    return highestOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillConnectionPathwayEnergyVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    const CaloHitList connectionPathwayU(LArConnectionPathwayHelper::GetConnectionPathwayProjection(pAlgorithm, protoShowerU, middleShowerStart3D, TPC_VIEW_U));
    const CaloHitList connectionPathwayV(LArConnectionPathwayHelper::GetConnectionPathwayProjection(pAlgorithm, protoShowerV, middleShowerStart3D, TPC_VIEW_V));
    const CaloHitList connectionPathwayW(LArConnectionPathwayHelper::GetConnectionPathwayProjection(pAlgorithm, protoShowerW, middleShowerStart3D, TPC_VIEW_W));

    electronTreeVariables.m_pathwayEnergyMeanU = LArConnectionPathwayHelper::CharacteriseEnergyMean(connectionPathwayU);
    electronTreeVariables.m_pathwayEnergyMeanV = LArConnectionPathwayHelper::CharacteriseEnergyMean(connectionPathwayV);
    electronTreeVariables.m_pathwayEnergyMeanW = LArConnectionPathwayHelper::CharacteriseEnergyMean(connectionPathwayW);

    electronTreeVariables.m_pathwayEnergySigmaU = LArConnectionPathwayHelper::CharacteriseEnergySigma(connectionPathwayU, electronTreeVariables.m_pathwayEnergyMeanU);
    electronTreeVariables.m_pathwayEnergySigmaV = LArConnectionPathwayHelper::CharacteriseEnergySigma(connectionPathwayV, electronTreeVariables.m_pathwayEnergyMeanV);
    electronTreeVariables.m_pathwayEnergySigmaW = LArConnectionPathwayHelper::CharacteriseEnergySigma(connectionPathwayW, electronTreeVariables.m_pathwayEnergyMeanW);

    electronTreeVariables.m_pathwayMinEnergyMean = std::min(std::min(electronTreeVariables.m_pathwayEnergyMeanU, electronTreeVariables.m_pathwayEnergyMeanV), 
        electronTreeVariables.m_pathwayEnergyMeanW);

    electronTreeVariables.m_pathwayMaxEnergyMean = std::max(std::max(electronTreeVariables.m_pathwayEnergyMeanU, electronTreeVariables.m_pathwayEnergyMeanV), 
        electronTreeVariables.m_pathwayEnergyMeanW);

    for (const float mean : {electronTreeVariables.m_pathwayEnergyMeanU, electronTreeVariables.m_pathwayEnergyMeanV, electronTreeVariables.m_pathwayEnergyMeanW})
    {
        if ((std::fabs(mean - electronTreeVariables.m_pathwayMinEnergyMean) > std::numeric_limits<float>::epsilon()) &&
            (std::fabs(mean - electronTreeVariables.m_pathwayMaxEnergyMean) > std::numeric_limits<float>::epsilon()))
        {
            electronTreeVariables.m_pathwayMiddleEnergyMean = mean;
        }
    }

    electronTreeVariables.m_pathwayMinEnergySigma = std::min(std::min(electronTreeVariables.m_pathwayEnergySigmaU, electronTreeVariables.m_pathwayEnergySigmaV), 
        electronTreeVariables.m_pathwayEnergySigmaW);

    electronTreeVariables.m_pathwayMaxEnergySigma = std::max(std::max(electronTreeVariables.m_pathwayEnergySigmaU, electronTreeVariables.m_pathwayEnergySigmaV), 
        electronTreeVariables.m_pathwayEnergySigmaW);

    for (const float sigma : {electronTreeVariables.m_pathwayEnergySigmaU, electronTreeVariables.m_pathwayEnergySigmaV, electronTreeVariables.m_pathwayEnergySigmaW})
    {
        if ((std::fabs(sigma - electronTreeVariables.m_pathwayMinEnergySigma) > std::numeric_limits<float>::epsilon()) &&
            (std::fabs(sigma - electronTreeVariables.m_pathwayMaxEnergySigma) > std::numeric_limits<float>::epsilon()))
        {
            electronTreeVariables.m_pathwayMiddleEnergySigma = sigma;
        }
    }

    electronTreeVariables.m_pathwayMaxEnergySigma = std::min(electronTreeVariables.m_pathwayMaxEnergySigma, 4.f);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArConnectionPathwayHelper::GetConnectionPathwayProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const CartesianVector &showerStart3D, const HitType hitType)
{
    const CartesianVector projectedShowerStart(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), showerStart3D, hitType));

    CartesianPointVector spinePositions;

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
        spinePositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult spineFit(&spinePositions, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

    CaloHitList projectedPathwayHitList;
    float lShowerStart(0.f), tShowerStart(0.f);
    spineFit.GetLocalPosition(projectedShowerStart, lShowerStart, tShowerStart);

    const bool isDownstream(protoShower.m_connectionPathway.m_startDirection.GetZ() > 0.f);

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        float thisL(0.f), thisT(0.f);
        spineFit.GetLocalPosition(hitPosition, thisL, thisT);

        if (isDownstream && (thisL < lShowerStart))
            projectedPathwayHitList.push_back(pCaloHit);

        if (!isDownstream && (thisL > lShowerStart))
            projectedPathwayHitList.push_back(pCaloHit);
    }

    return projectedPathwayHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::CharacteriseEnergyMean(const CaloHitList &caloHitList)
{
    float totalEnergy(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
        totalEnergy += (pCaloHit->GetElectromagneticEnergy() * 1000.f);

    return (caloHitList.size() == 0 ? -10.f : totalEnergy / static_cast<float>(caloHitList.size()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::CharacteriseEnergySigma(const CaloHitList &caloHitList, const float energyMean)
{ 
    float energySigma(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
        energySigma += std::pow((pCaloHit->GetElectromagneticEnergy() * 1000.f) - energyMean, 2);

    return (caloHitList.size() == 0 ? -10.f : std::sqrt(energySigma / static_cast<float>(caloHitList.size())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStarts3D(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerMatch &protoShowerMatch, const CartesianVector &nuVertexPosition, CartesianPointVector &showerStarts3D)
{
    const ElectronProtoShower &protoShowerU(protoShowerMatch.m_protoShowerU);
    const ElectronProtoShower &protoShowerV(protoShowerMatch.m_protoShowerV);
    const ElectronProtoShower &protoShowerW(protoShowerMatch.m_protoShowerW);
    const Consistency consistency(protoShowerMatch.m_consistencyType);

    // ISOBEL YOU ARE SO LAZY
    float maxSeparation(5.f);

    //std::cout << "consistency: " << consistency << std::endl;

    bool uFound(false), vFound(false), wFound(false);
    CartesianVector uShowerStart3D(0.f, 0.f, 0.f), vShowerStart3D(0.f, 0.f, 0.f), wShowerStart3D(0.f, 0.f, 0.f);

    CaloHitList caloHitList3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHitList3D);

    if (consistency == Consistency::POSITION)
    {
        LArConnectionPathwayHelper::FindShowerVertexFromPosition(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D);
        vShowerStart3D = uShowerStart3D;
        wShowerStart3D = uShowerStart3D;

        if (LArClusterHelper::GetClosestDistance(uShowerStart3D, caloHitList3D) < 3.f)
        {
            //std::cout << "FOUND FROM POSITION" << std::endl;
            uFound = true; vFound = true; wFound = true;
        }
    }
    else if (consistency == Consistency::DIRECTION)
    {
        if (LArConnectionPathwayHelper::FindShowerVertexFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, vShowerStart3D, wShowerStart3D))
        {
            if (LArClusterHelper::GetClosestDistance(uShowerStart3D, caloHitList3D) < 3.f)
            {
                //std::cout << "FOUND U FROM DIRECTION" << std::endl;
                uFound = true;
            }

            if (LArClusterHelper::GetClosestDistance(vShowerStart3D, caloHitList3D) < 3.f)
            {
                //std::cout << "FOUND V FROM DIRECTION" << std::endl;
                vFound = true;
            }

            if (LArClusterHelper::GetClosestDistance(wShowerStart3D, caloHitList3D) < 3.f)
            {
                //std::cout << "FOUND W FROM DIRECTION" << std::endl;
                wFound = true;
            }
        }
    }

    if (!uFound || (consistency == Consistency::X_PROJECTION))
    {
        if (LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxSeparation, uShowerStart3D) || 
            LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, pShowerPfo, protoShowerU, protoShowerV, protoShowerW, maxSeparation, uShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, pShowerPfo, protoShowerU, protoShowerW, protoShowerV, maxSeparation, uShowerStart3D))
        {
                uFound = true;
            /*
            if (LArClusterHelper::GetClosestDistance(uShowerStart3D, caloHitList3D) < 3.f)
            {
                //std::cout << "FOUND U FROM X PROJECTION" << std::endl;
                uFound = true;
            }
            */
        }

        if (LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxSeparation, vShowerStart3D) || 
            LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, pShowerPfo, protoShowerV, protoShowerW, protoShowerU, maxSeparation, vShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, pShowerPfo, protoShowerV, protoShowerU, protoShowerW, maxSeparation, vShowerStart3D))
        {
                vFound = true;
            /*
            if (LArClusterHelper::GetClosestDistance(vShowerStart3D, caloHitList3D) < 3.f)
            {
                //std::cout << "FOUND V FROM X PROJECTION" << std::endl;
                vFound = true;
            }
            */
        }

        if (LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, pShowerPfo, protoShowerW, protoShowerU, protoShowerV, maxSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, pShowerPfo, protoShowerW, protoShowerV, protoShowerU, maxSeparation, wShowerStart3D))
        {
                wFound = true;
            /*
            if (LArClusterHelper::GetClosestDistance(wShowerStart3D, caloHitList3D) < 3.f)
            {
                //std::cout << "FOUND W FROM X PROJECTION" << std::endl;
                wFound = true;
            }
            */
        }
    }

    CartesianPointVector tempShowerStarts3D;

    if (uFound)
        tempShowerStarts3D.push_back(uShowerStart3D);

    if (vFound)
        tempShowerStarts3D.push_back(vShowerStart3D);

    if (wFound)
        tempShowerStarts3D.push_back(wShowerStart3D);

    std::sort(tempShowerStarts3D.begin(), tempShowerStarts3D.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertexPosition));

    if (tempShowerStarts3D.empty())
        return false;

    showerStarts3D.push_back(tempShowerStarts3D.front());
    showerStarts3D.push_back((tempShowerStarts3D.size() == 3) ? tempShowerStarts3D[1] : tempShowerStarts3D[0]);
    showerStarts3D.push_back(tempShowerStarts3D.back());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Determining Whether ProtoShowers Are Consistent
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation)
{
    float metric(0.0);

    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    return LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, showerStartU, showerStartV, showerStartW, maxXSeparation, 
        maxSeparation, metric);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const CartesianVector &showerStartU, 
    const CartesianVector &showerStartV, const CartesianVector &showerStartW, const float maxXSeparation, const float maxSeparation, float &metric)
{
    const float xSeparationUV(std::fabs(showerStartU.GetX() - showerStartV.GetX()));
    const float xSeparationUW(std::fabs(showerStartU.GetX() - showerStartW.GetX()));
    const float xSeparationVW(std::fabs(showerStartV.GetX() - showerStartW.GetX()));

    if ((xSeparationUV > maxXSeparation) || (xSeparationUW > maxXSeparation) || (xSeparationVW > maxXSeparation))
        return false;

    float chi2(0.f);
    CartesianVector projectionU(0.f, 0.f, 0.f), projectionV(0.f, 0.f, 0.f), projectionW(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, showerStartV, showerStartW, projectionU, chi2);
    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, showerStartW, showerStartU, projectionV, chi2);
    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, showerStartU, showerStartV, projectionW, chi2);

    const float separationU((projectionU - showerStartU).GetMagnitude());
    const float separationV((projectionV - showerStartV).GetMagnitude());
    const float separationW((projectionW - showerStartW).GetMagnitude());

    /////////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionU, "U PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionV, "V PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionW, "W PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartU, "U SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartV, "V SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartW, "W SHOWER START", BLACK, 2);
    std::cout << "separationU: " << separationU << std::endl;
    std::cout << "separationV: " << separationV << std::endl;
    std::cout << "separationW: " << separationW << std::endl;
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    metric = (separationU + separationV + separationW) / 3.f;

    return (metric < maxSeparation);
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle)
{
    float metric(0.0);

    return LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxOpeningAngle, metric);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle, float &metric)
{
    const CartesianVector &directionU1(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &directionV1(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &directionW1(protoShowerW.m_connectionPathway.m_startDirection);

    if (LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU1, directionV1, directionW1, maxOpeningAngle, metric))
    {
        return true;
    }
    else
    {
        const bool isDownstream(protoShowerW.m_showerCore.m_startPosition.GetZ() > protoShowerW.m_connectionPathway.m_startPosition.GetZ());

        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
            spinePositionsU.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
            spinePositionsV.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
            spinePositionsW.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const CartesianVector directionU2(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionV2(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionW2(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        return LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU2, directionV2, directionW2, maxOpeningAngle, metric);
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// not by reference since we might change them...
bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, CartesianVector directionU, 
    CartesianVector directionV, CartesianVector directionW, const float maxOpeningAngle, float &metric)
{
    ////std::cout << "maxOpeningAngle: " << maxOpeningAngle << std::endl;

    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    float radians((2.f * M_PI) / 180.f);
    bool isIsochronous((wireDeviationU < radians) && (wireDeviationV < radians) && (wireDeviationW < radians));

    if (isIsochronous)
    {
        int positiveCount(0);
        positiveCount += directionU.GetX() > 0.f ? 0 : 1;
        positiveCount += directionV.GetX() > 0.f ? 0 : 1;
        positiveCount += directionW.GetX() > 0.f ? 0 : 1;

        if (positiveCount >= 2)
        {
            directionU = CartesianVector(std::fabs(directionU.GetX()), 0.f, directionU.GetZ());
            directionV = CartesianVector(std::fabs(directionV.GetX()), 0.f, directionV.GetZ());
            directionW = CartesianVector(std::fabs(directionW.GetX()), 0.f, directionW.GetZ());
        }
    }

    if (directionU.GetX() * directionV.GetX() < 0.f)
        return false;

    if (directionU.GetX() * directionW.GetX() < 0.f)
        return false;

    if (directionV.GetX() * directionW.GetX() < 0.f)
        return false;
    
    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    if (isIsochronous)
    {
        openingAngleU = std::min(openingAngleU, 180.f - openingAngleU);
        openingAngleV = std::min(openingAngleV, 180.f - openingAngleV);
        openingAngleW = std::min(openingAngleW, 180.f - openingAngleW);
    }

    /////////////////////////////////
    /*
    const CartesianVector &nuVertexU(protoShowerU.m_connectionPathway.m_startPosition);;
    const CartesianVector &nuVertexV(protoShowerV.m_connectionPathway.m_startPosition);
    const CartesianVector &nuVertexW(protoShowerW.m_connectionPathway.m_startPosition);
    const CartesianVector endU(nuVertexU + (directionU * 100.f));
    const CartesianVector endV(nuVertexV + (directionV * 100.f));
    const CartesianVector endW(nuVertexW + (directionW * 100.f));
    const CartesianVector projectionEndU(nuVertexU + (projectionU * 100.f));
    const CartesianVector projectionEndV(nuVertexV + (projectionV * 100.f));
    const CartesianVector projectionEndW(nuVertexW + (projectionW * 100.f));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexU, &projectionEndU, "PROJECTION U", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexV, &projectionEndV, "PROJECTION V", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexW, &projectionEndW, "PROJECTION W", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexU, &endU, "DIRECTION U", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexV, &endV, "DIRECTION V", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexW, &endW, "DIRECTION W", BLACK, 2, 2);
    std::cout << "angularDeviationU: " << openingAngleU << std::endl;
    std::cout << "angularDeviationV: " << openingAngleV << std::endl;
    std::cout << "angularDeviationW: " << openingAngleW << std::endl;
   //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
   */
    /////////////////////////////////

    if ((openingAngleU > maxOpeningAngle) || (openingAngleV > maxOpeningAngle) || (openingAngleW > maxOpeningAngle))
    {
        //std::cout << "OPENING ANGLE IS TOO BIG" << std::endl;
        return false;
    }

    metric = (openingAngleU + openingAngleV + openingAngleW) / 3.f;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Finding a 3D Direction
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindDirection3D(pandora::Algorithm *const pAlgorithm, const bool isDownstream, const CartesianVector &startPositionU, 
    const CartesianVector &startPositionV, const CartesianVector &startPositionW, const CartesianVector &startPosition3D, const CartesianVector &directionU, 
    const CartesianVector &directionV, const CartesianVector &directionW, CartesianVector &direction3D)
{
    float metric(0.f);
    if (!LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU, directionV, directionW, 5.f, metric))
    {
        //std::cout << "POST SHOWER START DIRECTIONS ARE NOT CONSISTENT" << std::endl;
        return false;
    }

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float uDisplacement(std::fabs(directionU.GetCosOpeningAngle(xAxis)));
    const float vDisplacement(std::fabs(directionV.GetCosOpeningAngle(xAxis)));
    const float wDisplacement(std::fabs(directionW.GetCosOpeningAngle(xAxis)));

    const CartesianVector directionSeedU(startPositionU + (directionU * uDisplacement));
    const CartesianVector directionSeedV(startPositionV + (directionV * vDisplacement));
    const CartesianVector directionSeedW(startPositionW + (directionW * wDisplacement));

    if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, directionSeedU, directionSeedV, directionSeedW, 5.f, 2.f, metric))
    {
        //std::cout << "POST SHOWER DIRECTION SEEDS ARE NOT CONSISTENT" << std::endl;
        return false;
    }

    CartesianVector directionSeed3D(0.f, 0.f, 0.f);

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, directionSeedU, directionSeedV, directionSeedW, directionSeed3D, chi2);

    direction3D = (directionSeed3D - startPosition3D);

    if (isDownstream && (direction3D.GetZ() < 0.f))
        direction3D = direction3D * (-1.f);

    if (!isDownstream && (direction3D.GetZ() > 0.f))
        direction3D = direction3D * (-1.f);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Finding 3D Shower Start Position
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromPosition(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    const CartesianVector showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector showerStartW(protoShowerW.m_showerCore.m_startPosition);

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV, showerStartW, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm, const CartesianVector nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    CartesianVector uShowerStart3D(0.f, 0.f, 0.f), vShowerStart3D(0.f, 0.f, 0.f), wShowerStart3D(0.f, 0.f, 0.f);

    if(!LArConnectionPathwayHelper::FindShowerVertexFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, 
        vShowerStart3D, wShowerStart3D))
    {
        return false;
    }

    float minSeparation(std::numeric_limits<float>::max());

    for (const CartesianVector showerStart : {uShowerStart3D, vShowerStart3D, wShowerStart3D})
    {
        const float separation((showerStart - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            minSeparation = separation;
            showerStart3D = showerStart;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm,
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &uShowerStart3D, 
    CartesianVector &vShowerStart3D, CartesianVector &wShowerStart3D)
{
    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    const CartesianVector &projectedNuVertexU(protoShowerU.m_connectionPathway.m_startPosition);
    const CartesianVector &projectedNuVertexV(protoShowerV.m_connectionPathway.m_startPosition);
    const CartesianVector &projectedNuVertexW(protoShowerW.m_connectionPathway.m_startPosition);

    const CartesianVector &peakDirectionU(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionV(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionW(protoShowerW.m_connectionPathway.m_startDirection);

    const CartesianVector &displacementU(showerStartU - projectedNuVertexU);
    const CartesianVector &displacementV(showerStartV - projectedNuVertexV);
    const CartesianVector &displacementW(showerStartW - projectedNuVertexW);

    const float transverseU(peakDirectionU.GetCrossProduct(displacementU).GetMagnitude());
    const float transverseV(peakDirectionV.GetCrossProduct(displacementV).GetMagnitude());
    const float transverseW(peakDirectionW.GetCrossProduct(displacementW).GetMagnitude());

    if ((transverseU > 1.f) || (transverseV > 1.f) || (transverseW > 1.f))
        return false;

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float cosThetaU(std::fabs(peakDirectionU.GetCosOpeningAngle(xAxis)));
    const float cosThetaV(std::fabs(peakDirectionV.GetCosOpeningAngle(xAxis)));
    const float cosThetaW(std::fabs(peakDirectionW.GetCosOpeningAngle(xAxis)));

    const float x1(showerStartU.GetX());
    const CartesianVector showerStartV1(projectedNuVertexV + (peakDirectionV * (std::fabs(x1 - projectedNuVertexV.GetX()) * cosThetaV)));
    const CartesianVector showerStartW1(projectedNuVertexW + (peakDirectionW * (std::fabs(x1 - projectedNuVertexW.GetX()) * cosThetaW)));

    const float x2(showerStartV.GetX());
    const CartesianVector showerStartU2(projectedNuVertexU + (peakDirectionU * (std::fabs(x2 - projectedNuVertexU.GetX()) * cosThetaU)));
    const CartesianVector showerStartW2(projectedNuVertexW + (peakDirectionW * (std::fabs(x2 - projectedNuVertexW.GetX()) * cosThetaW)));

    const float x3(showerStartW.GetX());
    const CartesianVector showerStartU3(projectedNuVertexU + (peakDirectionU * (std::fabs(x3 - projectedNuVertexU.GetX()) * cosThetaU)));
    const CartesianVector showerStartV3(projectedNuVertexV + (peakDirectionV * (std::fabs(x3 - projectedNuVertexV.GetX()) * cosThetaV)));

    float chi2(0.0);

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV1, showerStartW1, uShowerStart3D, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU2, showerStartV, showerStartW2, vShowerStart3D, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU3, showerStartV3, showerStartW, wShowerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// this is for gamma alg usage
bool LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector nuVertexPosition, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    const float maxSeparation, CartesianVector &showerStart3D)
{
    CartesianVector showerStartU(0.f, 0.f, 0.f), showerStartV(0.f, 0.f, 0.f), showerStartW(0.f, 0.f, 0.f);

    const bool uFound(LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxSeparation, showerStartU));

    const bool vFound(LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxSeparation, showerStartV));

    const bool wFound(LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxSeparation, showerStartW));

    bool found(false);
    float minSeparation(std::numeric_limits<float>::max());

    if (uFound)
    {
        const float separation((showerStartU - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            minSeparation = separation;
            showerStart3D = showerStartU;
        }
    }

    if (vFound)
    {
        const float separation((showerStartV - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            minSeparation = separation;
            showerStart3D = showerStartV;
        }
    }

    if (wFound)
    {
        const float separation((showerStartW - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            minSeparation = separation;
            showerStart3D = showerStartW;
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, CartesianVector &showerStart3D)
{
    std::cout << "YO IZZIE WE COULD DO A THING HERE WE'RE WE TRY AND SEE IF HIGHER OR LOWER L HYPOTHESIS IS BETTER" << std::endl;

    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f);

    //std::cout << "HERE 1" << std::endl;
    const HitType hitType(protoShower.m_spineHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_spineHitList.front()->GetHitType());
    const HitType hitType2(protoShower2.m_spineHitList.front()->GetHitType());
    //std::cout << "HERE 2" << std::endl;

    bool found1(false);
    float lowestL1(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower1.m_spineHitList)
    {
        if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart.GetX()) > 0.5f)
            continue;

        float lVertex(protoShower1.m_connectionPathway.m_startDirection.GetDotProduct(pCaloHit->GetPositionVector() - protoShower1.m_connectionPathway.m_startPosition));

        if ((lVertex > 0.f) && (lVertex < lowestL1))
        {
            lowestL1 = lVertex;
            showerStart1 = pCaloHit->GetPositionVector();
            found1 = true;
        }
    }

    //std::cout << "AAAAA" << std::endl;
    if (!found1)
        return false;
    //std::cout << "BBBBB" << std::endl;

    bool found2(false);
    float lowestL2(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower2.m_spineHitList)
    {
        if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart.GetX()) > 0.5f)
            continue;

        float lVertex(protoShower2.m_connectionPathway.m_startDirection.GetDotProduct(pCaloHit->GetPositionVector() - protoShower2.m_connectionPathway.m_startPosition));

        if ((lVertex > 0.f) && (lVertex < lowestL2))
        {
            lowestL2 = lVertex;
            showerStart2 = pCaloHit->GetPositionVector();
            found2 = true;
        }
    }
    //std::cout << "CCCCCC" << std::endl;
    if (!found2)
        return false;
    //std::cout << "DDDDDD" << std::endl;
    // Now make sure that they agree...

    float chi2(0.f);
    CartesianVector projection(0.f, 0.f, 0.f), projection1(0.f, 0.f, 0.f), projection2(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType1, hitType2, showerStart1, showerStart2, projection, chi2);
    const float separationSquared((projection - showerStart).GetMagnitudeSquared());

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType2, showerStart, showerStart2, projection1, chi2);
    const float separationSquared1((projection1 - showerStart1).GetMagnitudeSquared());

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);
    const float separationSquared2((projection2 - showerStart2).GetMagnitudeSquared());

    //std::cout << "EEEEEEEEEEEEE" << std::endl;

    //////////////////
    /*
    std::cout << "SETTING THE NEW VERTEX" << std::endl;
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projection, "PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projection1, "1 PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projection2, "2 PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart, "SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart1, "1 SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart2, "2 SHOWER START", BLACK, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////
    
    if ((separationSquared > maxSeparation * maxSeparation) || (separationSquared1 > maxSeparation * maxSeparation) || 
        (separationSquared2 > maxSeparation * maxSeparation))
    {
        return false;
    }
    //std::cout << "FFFFFFFFFFFF" << std::endl;

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), hitType, hitType1, hitType2, showerStart, showerStart1, showerStart2, showerStart3D, chi2);

    //std::cout << "GGGGGGGGGGGGGGG" << std::endl;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShower &protoShower, const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.m_spineHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_spineHitList.front()->GetHitType());

    bool found1(false);
    float lowestL1(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower1.m_spineHitList)
    {
        if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart.GetX()) > 0.5f)
            continue;

        float lVertex(protoShower1.m_connectionPathway.m_startDirection.GetDotProduct(pCaloHit->GetPositionVector() - protoShower1.m_connectionPathway.m_startPosition));

        if ((lVertex > 0.f) && (lVertex < lowestL1))
        {
            lowestL1 = lVertex;
            showerStart1 = pCaloHit->GetPositionVector();
            found1 = true;
        }
    }

    if (!found1)
        return false;

    // Now make sure that they agree...
    const CaloHitList &caloHitList2(protoShower2.m_spineHitList);

    float chi2(0.f);
    CartesianVector projection2(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);
    const float separation(LArClusterHelper::GetClosestDistance(projection2, caloHitList2));

    //////////////////

    //////////////////
    /*
    std::cout << "SETTING THE NEW VERTEX" << std::endl;
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projection2, "PROJECTION", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart, "SHOWER START", VIOLET, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart1, "1 SHOWER START", VIOLET, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////

    
    if (separation > maxSeparation)
        return false;

    LArGeometryHelper::MergeTwoPositions3D(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::GetMinMiddleMax(const float value1, const float value2, const float value3, float &minValue, float &middleValue, 
    float &maxValue)
{
    minValue = std::min(std::min(value1, value2), value3);
    maxValue = std::max(std::max(value1, value2), value3);
    middleValue = minValue;

    for (const float value : {value1, value2, value3})
    {
        if ((std::fabs(value - minValue) < std::numeric_limits<float>::epsilon()) || 
            (std::fabs(value - maxValue) < std::numeric_limits<float>::epsilon()))
        {
            continue;
        }

        middleValue = value;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::GetMinMiddleMax(const double value1, const double value2, const double value3, double &minValue, double &middleValue, 
    double &maxValue)
{
    minValue = std::min(std::min(value1, value2), value3);
    maxValue = std::max(std::max(value1, value2), value3);
    middleValue = minValue;

    for (const double value : {value1, value2, value3})
    {
        if ((std::fabs(value - minValue) < std::numeric_limits<double>::epsilon()) || 
            (std::fabs(value - maxValue) < std::numeric_limits<double>::epsilon()))
        {
            continue;
        }

        middleValue = value;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Sorting Functions
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::SortByDistanceToPoint::operator()(const CartesianVector &lhs, const CartesianVector &rhs)
{
    return (m_referencePoint.GetDistanceSquared(lhs) < m_referencePoint.GetDistanceSquared(rhs));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::SortByDistanceToPoint::operator()(const CaloHit *const lhs, const CaloHit *const rhs)
{
    return (m_referencePoint.GetDistanceSquared(lhs->GetPositionVector()) < m_referencePoint.GetDistanceSquared(rhs->GetPositionVector()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
