/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

#include "TMVA/Reader.h"

using namespace pandora;

namespace lar_content
{

ShowerStartRefinementAlgorithm::ShowerStartRefinementAlgorithm() : 
    m_binSize(0.005),     
    m_electronFraction(0.3f),
    m_createTrainingTrees(true),
    m_minElectronCompleteness(0.33f),
    m_minElectronPurity(0.5f),
    m_minGammaCompleteness(0.33f),
    m_thresholdSignalGammaDisplacement(-3.f),
    m_electronTMVACut(-0.2987),
    m_TMVAReader("")
{
    m_TMVAReader.AddVariable("PathwayLengthMin", &m_TMVAElectronTreeVariables.m_pathwayLengthMin);
    m_TMVAReader.AddVariable("PathwayMaxScatteringAngle", &m_TMVAElectronTreeVariables.m_pathwayMaxScatteringAngle);
    m_TMVAReader.AddVariable("PathwayMiddleEnergySigma", &m_TMVAElectronTreeVariables.m_pathwayMiddleEnergySigma);
    m_TMVAReader.AddVariable("MaxNPostShowerStartHits", &m_TMVAElectronTreeVariables.m_maxNPostShowerStartHits);
    m_TMVAReader.AddVariable("MaxPostShowerStartScatterAngle", &m_TMVAElectronTreeVariables.m_maxPostShowerStartScatterAngle);
    m_TMVAReader.AddVariable("PostShowerStartOpeningAngle", &m_TMVAElectronTreeVariables.m_postShowerStartOpeningAngle);
    m_TMVAReader.AddVariable("PostShowerStartMeanTransverseAngle", &m_TMVAElectronTreeVariables.m_postShowerStartMeanTransverseAngle);
    m_TMVAReader.AddVariable("PostShowerStartEnergyWeightedMeanRadialDistance", &m_TMVAElectronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance);
    m_TMVAReader.AddVariable("PostShowerStartEstimatedMoliereRadius", &m_TMVAElectronTreeVariables.m_postShowerStartEstimatedMoliereRadius);
    m_TMVAReader.AddVariable("PostShowerStartInitialGapSize", &m_TMVAElectronTreeVariables.m_postShowerStartInitialGapSize);
    m_TMVAReader.AddVariable("InitialRegionDistanceToNuVertex", &m_TMVAElectronTreeVariables.m_initialRegionDistanceToNuVertex);
    m_TMVAReader.AddVariable("NViewsWithAmbiguousHits", &m_TMVAElectronTreeVariables.m_nViewsWithAmbiguousHits);
    m_TMVAReader.AddVariable("AmbiguousHitMinUnaccountedEnergy", &m_TMVAElectronTreeVariables.m_ambiguousHitMinUnaccountedEnergy);

    std::string weightFilePath = "/dune/app/users/imawby/selection/connectionPathwayBDT/showerStartMiddle/dataset_electronsel/weights/TMVAClassification_BDTG.weights.xml";
    m_TMVAReader.BookMVA("BDTG", weightFilePath);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerStartRefinementAlgorithm::~ShowerStartRefinementAlgorithm()
{
    if (m_createTrainingTrees)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ElectronSignalTree", "ConnectionPathwayTrees.root", "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "THE LIMIT DOES NOT EXIST" << std::endl;
        }

        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ElectronBackgroundTree", "ConnectionPathwayTrees.root", "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "THE LIMIT DOES NOT EXIST" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::Run()
{
    PfoVector pfoVector;
    this->FillPfoVector(pfoVector);

    CartesianVector nuVertexPosition(0.f, 0.f, 0.f);
    if (this->GetNeutrinoVertex(nuVertexPosition) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    this->FillGammaHitMap();
    this->FillElectronHitMap();

    //this->InitialiseElectronTrees();

    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListU", pCaloHitListU));

    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListV", pCaloHitListV));

    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitListW));

    // run tools
    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        for (ShowerStartRefinementBaseTool *const pShowerStartRefinementTool : m_algorithmToolVector)
        {
            if (std::find(m_deletedPfos.begin(), m_deletedPfos.end(), pPfo) != m_deletedPfos.end())
                continue;

            pShowerStartRefinementTool->Run(this, pPfo, nuVertexPosition, pCaloHitListU, pCaloHitListV, pCaloHitListW);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillTree(const std::string &treeName, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayLengthMin", electronTreeVariables.m_pathwayLengthMin));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayLengthMiddle", electronTreeVariables.m_pathwayLengthMiddle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayLengthMax", electronTreeVariables.m_pathwayLengthMax));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayShowerStartDelta", electronTreeVariables.m_pathwayShowerStartDelta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMaxScatteringAngle", electronTreeVariables.m_pathwayMaxScatteringAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MinShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_minShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MiddleShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_middleShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MaxShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayEnergyMeanU", electronTreeVariables.m_pathwayEnergyMeanU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayEnergyMeanV", electronTreeVariables.m_pathwayEnergyMeanV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayEnergyMeanW", electronTreeVariables.m_pathwayEnergyMeanW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayEnergySigmaU", electronTreeVariables.m_pathwayEnergySigmaU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayEnergySigmaV", electronTreeVariables.m_pathwayEnergySigmaV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayEnergySigmaW", electronTreeVariables.m_pathwayEnergySigmaW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMinEnergyMean", electronTreeVariables.m_pathwayMinEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMiddleEnergyMean", electronTreeVariables.m_pathwayMiddleEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMaxEnergyMean", electronTreeVariables.m_pathwayMaxEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMinEnergySigma", electronTreeVariables.m_pathwayMinEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMiddleEnergySigma", electronTreeVariables.m_pathwayMiddleEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PathwayMaxEnergySigma", electronTreeVariables.m_pathwayMaxEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartLength", electronTreeVariables.m_postShowerStartLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartNHits", electronTreeVariables.m_postShowerStartNHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartNHitsU", electronTreeVariables.m_postShowerStartNHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartNHitsV", electronTreeVariables.m_postShowerStartNHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartNHitsW", electronTreeVariables.m_postShowerStartNHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MinNPostShowerStartHits", electronTreeVariables.m_minNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MaxNPostShowerStartHits", electronTreeVariables.m_maxNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartScatterAngle", electronTreeVariables.m_postShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartScatterAngleU", electronTreeVariables.m_postShowerStartScatterAngleU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartScatterAngleV", electronTreeVariables.m_postShowerStartScatterAngleV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartScatterAngleW", electronTreeVariables.m_postShowerStartScatterAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MinPostShowerStartScatterAngle", electronTreeVariables.m_minPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MaxPostShowerStartScatterAngle", electronTreeVariables.m_maxPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartMinHalfOpeningAngle", electronTreeVariables.m_postShowerStartMinHalfOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartMaxHalfOpeningAngle", electronTreeVariables.m_postShowerStartMaxHalfOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartOpeningAngle", electronTreeVariables.m_postShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartOpeningAngleAsymmetry", electronTreeVariables.m_postShowerStartOpeningAngleAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartMeanTransverseAngle", electronTreeVariables.m_postShowerStartMeanTransverseAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartMeanLWeightedTransverseAngle", electronTreeVariables.m_postShowerStartMeanLWeightedTransverseAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartMeanRadialDistance", electronTreeVariables.m_postShowerStartMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartRadialDistanceSigma", electronTreeVariables.m_postShowerStartRadialDistanceSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartEnergyWeightedMeanRadialDistance", electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartEstimatedMoliereRadius", electronTreeVariables.m_postShowerStartEstimatedMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartLWeightedMeanRadialDistance", electronTreeVariables.m_postShowerStartLWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartLWeightedRadialDistanceSigma", electronTreeVariables.m_postShowerStartLWeightedRadialDistanceSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartInitialGapSize", electronTreeVariables.m_postShowerStartInitialGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "PostShowerStartMaxGapSize", electronTreeVariables.m_postShowerStartMaxGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "InitialRegionDistanceToNuVertex", electronTreeVariables.m_initialRegionDistanceToNuVertex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "InitialRegionDistanceInGaps", electronTreeVariables.m_initialRegionDistanceInGaps));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "InitialRegionMaxGapSize", electronTreeVariables.m_initialRegionMaxGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "NViewsWithAmbiguousHits", electronTreeVariables.m_nViewsWithAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "NAmbiguousHits2D", electronTreeVariables.m_nAmbiguousHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MinNAmbiguousHits", electronTreeVariables.m_minNAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "MaxNAmbiguousHits", electronTreeVariables.m_maxNAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitUnaccountedEnergyU", electronTreeVariables.m_ambiguousHitUnaccountedEnergyU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitUnaccountedEnergyV", electronTreeVariables.m_ambiguousHitUnaccountedEnergyV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitUnaccountedEnergyW", electronTreeVariables.m_ambiguousHitUnaccountedEnergyW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitMinUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitMaxUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitShowerEnergyRatioU", electronTreeVariables.m_ambiguousHitShowerEnergyRatioU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitShowerEnergyRatioV", electronTreeVariables.m_ambiguousHitShowerEnergyRatioV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitShowerEnergyRatioW", electronTreeVariables.m_ambiguousHitShowerEnergyRatioW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitMinShowerEnergyRatio", electronTreeVariables.m_ambiguousHitMinShowerEnergyRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronSignalTree", "AmbiguousHitMaxShowerEnergyRatio", electronTreeVariables.m_ambiguousHitMaxShowerEnergyRatio));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayLengthMin", electronTreeVariables.m_pathwayLengthMin));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayLengthMiddle", electronTreeVariables.m_pathwayLengthMiddle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayLengthMax", electronTreeVariables.m_pathwayLengthMax));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayShowerStartDelta", electronTreeVariables.m_pathwayShowerStartDelta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMaxScatteringAngle", electronTreeVariables.m_pathwayMaxScatteringAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MinShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_minShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MiddleShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_middleShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MaxShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayEnergyMeanU", electronTreeVariables.m_pathwayEnergyMeanU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayEnergyMeanV", electronTreeVariables.m_pathwayEnergyMeanV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayEnergyMeanW", electronTreeVariables.m_pathwayEnergyMeanW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayEnergySigmaU", electronTreeVariables.m_pathwayEnergySigmaU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayEnergySigmaV", electronTreeVariables.m_pathwayEnergySigmaV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayEnergySigmaW", electronTreeVariables.m_pathwayEnergySigmaW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMinEnergyMean", electronTreeVariables.m_pathwayMinEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMiddleEnergyMean", electronTreeVariables.m_pathwayMiddleEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMaxEnergyMean", electronTreeVariables.m_pathwayMaxEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMinEnergySigma", electronTreeVariables.m_pathwayMinEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMiddleEnergySigma", electronTreeVariables.m_pathwayMiddleEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PathwayMaxEnergySigma", electronTreeVariables.m_pathwayMaxEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartLength", electronTreeVariables.m_postShowerStartLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartNHits", electronTreeVariables.m_postShowerStartNHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartNHitsU", electronTreeVariables.m_postShowerStartNHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartNHitsV", electronTreeVariables.m_postShowerStartNHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartNHitsW", electronTreeVariables.m_postShowerStartNHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MinNPostShowerStartHits", electronTreeVariables.m_minNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MaxNPostShowerStartHits", electronTreeVariables.m_maxNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartScatterAngle", electronTreeVariables.m_postShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartScatterAngleU", electronTreeVariables.m_postShowerStartScatterAngleU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartScatterAngleV", electronTreeVariables.m_postShowerStartScatterAngleV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartScatterAngleW", electronTreeVariables.m_postShowerStartScatterAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MinPostShowerStartScatterAngle", electronTreeVariables.m_minPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MaxPostShowerStartScatterAngle", electronTreeVariables.m_maxPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartMinHalfOpeningAngle", electronTreeVariables.m_postShowerStartMinHalfOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartMaxHalfOpeningAngle", electronTreeVariables.m_postShowerStartMaxHalfOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartOpeningAngle", electronTreeVariables.m_postShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartOpeningAngleAsymmetry", electronTreeVariables.m_postShowerStartOpeningAngleAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartMeanTransverseAngle", electronTreeVariables.m_postShowerStartMeanTransverseAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartMeanLWeightedTransverseAngle", electronTreeVariables.m_postShowerStartMeanLWeightedTransverseAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartMeanRadialDistance", electronTreeVariables.m_postShowerStartMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartRadialDistanceSigma", electronTreeVariables.m_postShowerStartRadialDistanceSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartEnergyWeightedMeanRadialDistance", electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartEstimatedMoliereRadius", electronTreeVariables.m_postShowerStartEstimatedMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartLWeightedMeanRadialDistance", electronTreeVariables.m_postShowerStartLWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartLWeightedRadialDistanceSigma", electronTreeVariables.m_postShowerStartLWeightedRadialDistanceSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartInitialGapSize", electronTreeVariables.m_postShowerStartInitialGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "PostShowerStartMaxGapSize", electronTreeVariables.m_postShowerStartMaxGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "InitialRegionDistanceToNuVertex", electronTreeVariables.m_initialRegionDistanceToNuVertex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "InitialRegionDistanceInGaps", electronTreeVariables.m_initialRegionDistanceInGaps));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "InitialRegionMaxGapSize", electronTreeVariables.m_initialRegionMaxGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "NViewsWithAmbiguousHits", electronTreeVariables.m_nViewsWithAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "NAmbiguousHits2D", electronTreeVariables.m_nAmbiguousHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MinNAmbiguousHits", electronTreeVariables.m_minNAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "MaxNAmbiguousHits", electronTreeVariables.m_maxNAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitUnaccountedEnergyU", electronTreeVariables.m_ambiguousHitUnaccountedEnergyU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitUnaccountedEnergyV", electronTreeVariables.m_ambiguousHitUnaccountedEnergyV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitUnaccountedEnergyW", electronTreeVariables.m_ambiguousHitUnaccountedEnergyW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitMinUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitMaxUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitShowerEnergyRatioU", electronTreeVariables.m_ambiguousHitShowerEnergyRatioU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitShowerEnergyRatioV", electronTreeVariables.m_ambiguousHitShowerEnergyRatioV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitShowerEnergyRatioW", electronTreeVariables.m_ambiguousHitShowerEnergyRatioW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitMinShowerEnergyRatio", electronTreeVariables.m_ambiguousHitMinShowerEnergyRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ElectronBackgroundTree", "AmbiguousHitMaxShowerEnergyRatio", electronTreeVariables.m_ambiguousHitMaxShowerEnergyRatio));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillPfoVector(PfoVector &pfoVector)
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        if (!pPfoList || pPfoList->empty())
        {
            std::cout << "ShowerStartRefinementAlgorithm: unable to find pfo list " << pfoListName << std::endl;
            continue;
        }

        pfoVector.insert(pfoVector.begin(), pPfoList->begin(), pPfoList->end());
    }

    // This ordering is important.
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &neutrinoVertex)
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    neutrinoVertex = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList ShowerStartRefinementAlgorithm::GetAllHitsOfType(const HitType hitType)
{
    CaloHitList viewHitList;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return viewHitList;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (pCaloHit->GetHitType() == hitType)
            viewHitList.push_back(pCaloHit);
    }

    return viewHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList ShowerStartRefinementAlgorithm::GetXIntervalHitsOfType(const ParticleFlowObject *const pShowerPfo, const HitType hitType)
{
    CaloHitList intervalHitList;

    ClusterList clustersU, clustersV, clustersW;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_U, clustersU); 
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_V, clustersV); 
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_W, clustersW); 

    CartesianVector uMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    CartesianVector uMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    CartesianVector vMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    CartesianVector vMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    CartesianVector wMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    CartesianVector wMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    if (!clustersU.empty())
        LArClusterHelper::GetClusterBoundingBox(clustersU.front(), uMin, uMax);

    if (!clustersV.empty())
        LArClusterHelper::GetClusterBoundingBox(clustersV.front(), vMin, vMax);

    if (!clustersW.empty())
        LArClusterHelper::GetClusterBoundingBox(clustersW.front(), wMin, wMax);

    float xMin(std::min(std::min(uMin.GetX(), vMin.GetX()), wMin.GetX()));
    float xMax(std::max(std::max(uMax.GetX(), vMax.GetX()), wMax.GetX()));

    CaloHitList viewHitList;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return intervalHitList;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (pCaloHit->GetHitType() != hitType)
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if ((hitPosition.GetX() > xMin) && (hitPosition.GetX() < xMax))
            intervalHitList.push_back(pCaloHit);
    }

    return intervalHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillOwnershipMaps()
{
    m_hitToClusterMapU.clear(); m_hitToClusterMapV.clear(); m_hitToClusterMapW.clear();
    m_clusterToPfoMapU.clear(); m_clusterToPfoMapV.clear(); m_clusterToPfoMapW.clear();

    // First fill pfo maps
    PfoVector pfoVector;

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        //PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, );
        PandoraContentApi::GetList(*this, pfoListName, pPfoList);

        if (!pPfoList || pPfoList->empty())
        {
            std::cout << "ShowerStartRefinementAlgorithm: unable to find pfo list " << pfoListName << std::endl;
            continue;
        }

        pfoVector.insert(pfoVector.begin(), pPfoList->begin(), pPfoList->end());
    }


    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        ClusterList twoDClusterList;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, twoDClusterList);

        for (const Cluster *const pCluster : twoDClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

            for (const CaloHit *const pIsolated : isolated)
            {
                if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                    caloHitList.push_back(pIsolated);
            }

            const HitType hitType(caloHitList.front()->GetHitType());
            HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
            ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

            for (const CaloHit *const pCaloHit : caloHitList)
                hitToClusterMap[pCaloHit] = pCluster;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }

    // Now fill cluster maps
    StringVector clusterListNames;
    clusterListNames.push_back("ClustersU");
    clusterListNames.push_back("ClustersV");
    clusterListNames.push_back("ClustersW");

    ClusterList clusterList;
    for (const std::string &clusterListName : clusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            std::cout << "ShowerStartRefinementAlgorithm: No cluster list found, returning..." << std::endl;
            throw;
        }

        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
            continue;

        if (!pCluster->IsAvailable())
        {
            std::cout << "CLUSTER IS NOT AVAILABLE ISOBE" << std::endl;
            throw;
        }

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolated)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        for (const CaloHit *const pCaloHit : caloHitList)
            hitToClusterMap[pCaloHit] = pCluster;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillGammaHitMap()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            
            if (pMCParticle->GetParticleId() == 22)
                m_gammaHitMap[pMCParticle].push_back(pCaloHit);
        }
        catch (...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillElectronHitMap()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        MCParticleVector contributingMCParticleVector;
        const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());

        for (const auto &mapEntry : weightMap)
            contributingMCParticleVector.push_back(mapEntry.first);

        std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

        float highestWeight(0.f);
        const MCParticle *highestElectronContributor(nullptr);

        for (const MCParticle *const pMCParticle : contributingMCParticleVector)
        {
            const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) && (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

            if (isLeadingElectron)
            {
                const float weight(weightMap.at(pMCParticle));

                if (weight > highestWeight)
                {
                    highestWeight = weight;
                    highestElectronContributor = pMCParticle;
                }
            }
        }

        if (highestElectronContributor)
            m_electronHitMap[highestElectronContributor].push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementAlgorithm::IsElectron(const ParticleFlowObject *const pPfo) const
{
    MCParticleVector mcElectronVector;

    for (auto &entry : m_electronHitMap)
        mcElectronVector.push_back(entry.first);

    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

    for (const MCParticle *const pMCElectron : mcElectronVector)
    {
        const CaloHitList &mcElectronHitList(m_electronHitMap.at(pMCElectron));
        const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcElectronHitList));

        const float completeness(static_cast<float>(sharedHitList.size()) / static_cast<float>(mcElectronHitList.size()));
        const float purity(static_cast<float>(sharedHitList.size()) / static_cast<float>(pfoHitList.size()));

        if (completeness < m_minElectronCompleteness)
            continue;

        if (purity < m_minElectronPurity)
            continue;

        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementAlgorithm::TMVAIsElectron(LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    m_TMVAElectronTreeVariables.m_pathwayLengthMin = electronTreeVariables.m_pathwayLengthMin;
    m_TMVAElectronTreeVariables.m_pathwayMaxScatteringAngle = electronTreeVariables.m_pathwayMaxScatteringAngle;
    m_TMVAElectronTreeVariables.m_pathwayMiddleEnergySigma = electronTreeVariables.m_pathwayMiddleEnergySigma;
    m_TMVAElectronTreeVariables.m_maxNPostShowerStartHits = electronTreeVariables.m_maxNPostShowerStartHits;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartScatterAngle = electronTreeVariables.m_maxPostShowerStartScatterAngle;
    m_TMVAElectronTreeVariables.m_postShowerStartOpeningAngle = electronTreeVariables.m_postShowerStartOpeningAngle;
    m_TMVAElectronTreeVariables.m_postShowerStartMeanTransverseAngle = electronTreeVariables.m_postShowerStartMeanTransverseAngle;
    m_TMVAElectronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance = electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance;
    m_TMVAElectronTreeVariables.m_postShowerStartEstimatedMoliereRadius = electronTreeVariables.m_postShowerStartEstimatedMoliereRadius;
    m_TMVAElectronTreeVariables.m_postShowerStartInitialGapSize = electronTreeVariables.m_postShowerStartInitialGapSize;
    m_TMVAElectronTreeVariables.m_initialRegionDistanceToNuVertex = electronTreeVariables.m_initialRegionDistanceToNuVertex;
    m_TMVAElectronTreeVariables.m_nViewsWithAmbiguousHits = electronTreeVariables.m_nViewsWithAmbiguousHits;
    m_TMVAElectronTreeVariables.m_ambiguousHitMinUnaccountedEnergy = electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy;

    float bdtScore(m_TMVAReader.EvaluateMVA("BDTG"));

    std::cout << "bdtScore: " << bdtScore << std::endl;

    return bdtScore > m_electronTMVACut;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// our signal here is actually gammas that have made a mistake by getting back to the nu vertex. i.e. their reco vertex is closer to the nu vertex than the truth says

bool ShowerStartRefinementAlgorithm::IsGamma(const ParticleFlowObject *const pPfo, const CartesianVector &nuVertexPosition) const
{
    const MCParticle *pMainMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    int pdg(std::abs(pMainMCParticle->GetParticleId()));

    if (pdg == 11)
        return false;

    std::cout << "DDDDDDDDD" << std::endl;

    // but does it contain a large chunk of a gamma even if it is not the main owner?
    // find gamma with highest number of shared hits that has a completeness of >50 % in one view
    if (pdg != 22)
    {
        unsigned int highestSharedHits(0);

        MCParticleVector mcGammaVector;
        for (auto &entry : m_gammaHitMap)
        {
            if (entry.first->GetParticleId() == 22)
            {
                const CaloHitList &mcGammaHitList(m_gammaHitMap.at(entry.first));

                if (mcGammaHitList.size() > 100)
                    mcGammaVector.push_back(entry.first);
            }
        }

        std::sort(mcGammaVector.begin(), mcGammaVector.end(), LArMCParticleHelper::SortByMomentum);

        CaloHitList pfoHitList;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

        for (const MCParticle *const pMCGamma : mcGammaVector)
        {
            const CaloHitList &mcGammaHitList(m_gammaHitMap.at(pMCGamma));
            const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcGammaHitList));

           // get each view completeness
            for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                const float completeness(static_cast<float>(LArMonitoringHelper::CountHitsByType(hitType, sharedHitList)) /
                    static_cast<float>(LArMonitoringHelper::CountHitsByType(hitType, mcGammaHitList)));

                if ((completeness > m_minGammaCompleteness) && (sharedHitList.size() > highestSharedHits))
                {
                    pMainMCParticle = pMCGamma;
                    highestSharedHits = sharedHitList.size();
                }
            }
        }
    }

    pdg = std::abs(pMainMCParticle->GetParticleId());

    if (pdg != 22)
        return false;

    std::cout << "AAAAAAAAAAAAA" << std::endl;

    CaloHitList caloHitList3D;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList3D);

    if (caloHitList3D.empty())
        return false;

    std::cout << "BBBBBBBBBBBBB" << std::endl;

    float closestDistance(std::numeric_limits<float>::max());
    CartesianVector showerVertex(0.f, 0.f, 0.f);

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            float jam = (pCaloHit->GetPositionVector() - nuVertexPosition).GetMagnitude();
            if (jam < closestDistance)
            {
                closestDistance = jam;
                showerVertex = pCaloHit->GetPositionVector();
            }
        }


        /*
        LArPointingCluster pointingCluster(clusters3D.front(), 20, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        const CartesianVector innerVertex(pointingCluster.GetInnerVertex().GetPosition());
        const CartesianVector outerVertex(pointingCluster.GetOuterVertex().GetPosition());

        if ((outerVertex - nuVertexPosition).GetMagnitude() < (innerVertex - nuVertexPosition).GetMagnitude())
        {
            showerVertex = outerVertex;
        }
        else
        {
            showerVertex = innerVertex;
        }
        */

    const CartesianVector &mcVertex(pMainMCParticle->GetVertex());
    const float mcSeparation((mcVertex - nuVertexPosition).GetMagnitude());
    const float recoSeparation((showerVertex - nuVertexPosition).GetMagnitude());
    const float difference(recoSeparation - mcSeparation);


    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertex, "mcVertex", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerVertex, "showerVertex", GREEN, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    std::cout << "CCCCCCCC" << std::endl;

    if (difference > m_thresholdSignalGammaDisplacement)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------


void ShowerStartRefinementAlgorithm::AddElectronPathway(const ParticleFlowObject *const pShowerPfo, const CaloHitList &pathwayHitList)
{
    // This is so incredibly lazy isobel
    this->FillOwnershipMaps();

    const HitType hitType(pathwayHitList.front()->GetHitType());

    ClusterList showerClusters2D, showerClusters3D;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, showerClusters2D);
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, showerClusters3D);

    const ThreeDSlidingFitResult showerSlidingFit(showerClusters3D.front(), 1000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const CartesianVector &showerDirection3D(showerSlidingFit.GetGlobalMaxLayerDirection());

    std::map<const ParticleFlowObject*, int> showerHitCountMap;

    for (const CaloHit *const pPathwayHit : pathwayHitList)
    {
        const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (hitToClusterMap.find(pPathwayHit) == hitToClusterMap.end())
            continue;

        if (clusterToPfoMap.find(hitToClusterMap.at(pPathwayHit)) == clusterToPfoMap.end())
            continue;

        const ParticleFlowObject *const pPathwayShower(clusterToPfoMap.at(hitToClusterMap.at(pPathwayHit)));

        if (LArPfoHelper::IsTrack(pPathwayShower))
            continue;

        if (pPathwayShower == pShowerPfo)
            continue;

        if (showerHitCountMap.find(pPathwayShower) == showerHitCountMap.end())
            showerHitCountMap[pPathwayShower] = 1;
        else
            showerHitCountMap[pPathwayShower] = showerHitCountMap[pPathwayShower] + 1;
    }

    PfoList significantShowersToMerge;

    for (const auto &entry : showerHitCountMap)
    {
        float contaminationRatio(static_cast<float>(entry.second) / static_cast<float>(pathwayHitList.size()));

        if (contaminationRatio < 0.3f)
            continue;

        ClusterList pathwayClusters3D;
        LArPfoHelper::GetClusters(entry.first, TPC_3D, pathwayClusters3D);

        if (pathwayClusters3D.size() == 0)
            continue;

        try
        {
            const ThreeDSlidingFitResult pathwaySlidingFit(pathwayClusters3D.front(), 1000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const CartesianVector &pathwayDirection(pathwaySlidingFit.GetGlobalMaxLayerDirection());

            const float openingAngle(pathwayDirection.GetOpeningAngle(showerDirection3D) * 180 / M_PI);

            if (openingAngle < 5.f)
                significantShowersToMerge.push_back(entry.first);
        }
        catch (...)
        {
        }
    }

    // Add in hits first, then deal with merges
    for (const CaloHit *const pPathwayHit : pathwayHitList)
    {
        const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);
        std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        const Cluster *pParentCluster(nullptr);
        const ParticleFlowObject *pParentPfo(nullptr);

        if (hitToClusterMap.find(pPathwayHit) != hitToClusterMap.end())
        {
            pParentCluster = hitToClusterMap.at(pPathwayHit);

            if (clusterToPfoMap.find(pParentCluster) != clusterToPfoMap.end())
            {
                pParentPfo = clusterToPfoMap.at(pParentCluster);

                if (pParentPfo == pShowerPfo)
                    continue;

                if (std::find(significantShowersToMerge.begin(), significantShowersToMerge.end(), pParentPfo) != significantShowersToMerge.end())
                    continue;
            }
        }

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

            CaloHitList clusterNormalHitList; const CaloHitList clusterIsolatedHitList(pParentCluster->GetIsolatedCaloHitList());
            pParentCluster->GetOrderedCaloHitList().FillCaloHitList(clusterNormalHitList);

            const bool isIsolated(std::find(clusterIsolatedHitList.begin(), clusterIsolatedHitList.end(), pPathwayHit) != clusterIsolatedHitList.end());

            if (!isIsolated && (clusterNormalHitList.size() == 1) && !(clusterIsolatedHitList.empty()))
            {
                const HitType isolatedHitType(LArClusterHelper::GetClusterHitType(pParentCluster));
                HitToClusterMap &isolatedHitToClusterMap(isolatedHitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);

                for (const CaloHit * const pIsolatedHit : clusterIsolatedHitList)
                {
                    isolatedHitToClusterMap.erase(pIsolatedHit);
                    const StatusCode isolatedStatusCode(PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pIsolatedHit));

                    if (isolatedStatusCode != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ISOBEL CANNOT REMOVE ISOLATED HIT?" << std::endl;
                        throw;
                    }
                }
            }

            const StatusCode statusCodeCluster(isIsolated ? PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pPathwayHit) : 
                PandoraContentApi::RemoveFromCluster(*this, pParentCluster, pPathwayHit));

            if (statusCodeCluster != STATUS_CODE_SUCCESS)
            {
                if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                {
                    std::cout << "ElectronStartRefinementTool: cluster jam" << std::endl;
                    throw StatusCodeException(statusCodeCluster);
                }

                if (pParentPfo)
                {
                    const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pParentPfo, pParentCluster));
                    const unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pParentPfo));

                    if (nHits == 0)
                        std::cout << "ElectronStartRefinementTool: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                    if (statusCodePfo != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ElectronStartRefinementTool: pfo jam" << std::endl;
                        throw;
                    }
                }

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pParentCluster));
            }
        }

        if (!PandoraContentApi::IsAvailable(*this, pPathwayHit))
        {
            std::cout << "CALO HIT IS NOT AVAILABLE!!" << std::endl;
            throw;
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, showerClusters2D.front(), pPathwayHit));
    }

    // Now handle the shower merges
    for (const ParticleFlowObject *const pShowerToMerge : significantShowersToMerge)
    {
        m_deletedPfos.push_back(pShowerToMerge);
        this->MergeAndDeletePfos(pShowerPfo, pShowerToMerge);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::SetElectronMetadata(const CartesianVector &nuVertexPosition, const ParticleFlowObject *const pShowerPfo)
{
    object_creation::ParticleFlowObject::Metadata metadata;
    metadata.m_propertiesToAdd["ShowerVertexX"] = nuVertexPosition.GetX();
    metadata.m_propertiesToAdd["ShowerVertexY"] = nuVertexPosition.GetY();
    metadata.m_propertiesToAdd["ShowerVertexZ"] = nuVertexPosition.GetZ();
    metadata.m_propertiesToAdd["dEdX"] = 2.3;
    metadata.m_propertiesToAdd["ActiveHybridElectronAlg"] = 1.f;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pShowerPfo, metadata));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::SetGammaVertex(const CartesianVector &showerVertex, const ParticleFlowObject *const pShowerPfo)
{
    if (!pShowerPfo->GetVertexList().empty())
    {
        if (pShowerPfo->GetVertexList().size() != 1)
        {
            std::cout << "vertex not equal to one!!" << std::endl;
            throw;
        }

        const Vertex *const pVertex(pShowerPfo->GetVertexList().front());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pVertex, "GammaVertices"));
    }


    const VertexList *pVertexList = NULL;
    std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = showerVertex;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, "GammaVertices"));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pShowerPfo, pVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::RemoveConnectionPathway(const ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower)
{
    const HitType hitType(protoShower.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, clusterList);

    std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

    for (const CaloHit *const pCaloHit : protoShower.m_connectionPathway.m_pathwayHitList)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, clusterList.front(), pCaloHit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ShowerStartRefinementTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        ShowerStartRefinementBaseTool *const pShowerStartRefinementTool(dynamic_cast<ShowerStartRefinementBaseTool *>(*iter));

        if (!pShowerStartRefinementTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pShowerStartRefinementTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "BinSize", m_binSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ElectronFraction", m_electronFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ThresholdSignalGammaDisplacement", m_thresholdSignalGammaDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ElectronTMVACut", m_electronTMVACut));

    PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

    /*
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerDistribution", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "THE LIMIT DOES NOT EXIST" << std::endl;
    }
    */
