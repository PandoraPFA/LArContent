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
    m_createTrainingTrees(false),
    m_hybridMode(true),
    m_electronTMVACut(-0.1),
    m_gammaTMVACut(-0.1),
    m_minElectronCompleteness(0.33f),
    m_minElectronPurity(0.5f),
    m_minGammaCompleteness(0.33f),
    m_thresholdSignalGammaDisplacement(-3.f),
    m_TMVAReader("")
{
    m_TMVAReader.AddVariable("PathwayLengthMin", &m_TMVAElectronTreeVariables.m_pathwayLengthMin);
    m_TMVAReader.AddVariable("MaxShowerStartPathwayScatteringAngle2D", &m_TMVAElectronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D);
    m_TMVAReader.AddVariable("MaxNPostShowerStartHits", &m_TMVAElectronTreeVariables.m_maxNPostShowerStartHits);
    m_TMVAReader.AddVariable("MaxPostShowerStartScatterAngle", &m_TMVAElectronTreeVariables.m_maxPostShowerStartScatterAngle);
    m_TMVAReader.AddVariable("MaxPostShowerStartNuVertexEnergyAsymmetry", &m_TMVAElectronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry);
    m_TMVAReader.AddVariable("MaxPostShowerStartShowerStartEnergyAsymmetry", &m_TMVAElectronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry);
    m_TMVAReader.AddVariable("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", &m_TMVAElectronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
    m_TMVAReader.AddVariable("MinPostShowerStartShowerStartMoliereRadius", &m_TMVAElectronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius);
    m_TMVAReader.AddVariable("MaxPostShowerStartOpeningAngle", &m_TMVAElectronTreeVariables.m_maxPostShowerStartOpeningAngle);
    m_TMVAReader.AddVariable("MaxFoundHitRatio", &m_TMVAElectronTreeVariables.m_maxFoundHitRatio);
    m_TMVAReader.AddVariable("MaxInitialGapSize", &m_TMVAElectronTreeVariables.m_maxInitialGapSize);
    m_TMVAReader.AddVariable("MinLargestProjectedGapSize", &m_TMVAElectronTreeVariables.m_minLargestProjectedGapSize);
    m_TMVAReader.AddVariable("NViewsWithAmbiguousHits", &m_TMVAElectronTreeVariables.m_nViewsWithAmbiguousHits);
    m_TMVAReader.AddVariable("AmbiguousHitMaxUnaccountedEnergy", &m_TMVAElectronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy);

    std::string weightFilePath = "TMVAClassification_BDTG_StandardVertex.weights.xml";
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

        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "GammaSignalTree", "ConnectionPathwayTrees.root", "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "THE LIMIT DOES NOT EXIST" << std::endl;
        }

        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "GammaBackgroundTree", "ConnectionPathwayTrees.root", "UPDATE"));
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
    std::cout << "electronTreeVariables.m_pathwayLengthMin: " << electronTreeVariables.m_pathwayLengthMin << std::endl;
    std::cout << "electronTreeVariables.m_initialRegionDistanceToNuVertex: " << electronTreeVariables.m_initialRegionDistanceToNuVertex << std::endl;
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "NConnectionPathways", electronTreeVariables.m_nConnectionPathways));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayLengthMin", electronTreeVariables.m_pathwayLengthMin));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayLengthMiddle", electronTreeVariables.m_pathwayLengthMiddle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayLengthMax", electronTreeVariables.m_pathwayLengthMax));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayShowerStartDelta", electronTreeVariables.m_pathwayShowerStartDelta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMaxScatteringAngle", electronTreeVariables.m_pathwayMaxScatteringAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_minShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_middleShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayEnergyMeanU", electronTreeVariables.m_pathwayEnergyMeanU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayEnergyMeanV", electronTreeVariables.m_pathwayEnergyMeanV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayEnergyMeanW", electronTreeVariables.m_pathwayEnergyMeanW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayEnergySigmaU", electronTreeVariables.m_pathwayEnergySigmaU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayEnergySigmaV", electronTreeVariables.m_pathwayEnergySigmaV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayEnergySigmaW", electronTreeVariables.m_pathwayEnergySigmaW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMinEnergyMean", electronTreeVariables.m_pathwayMinEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMiddleEnergyMean", electronTreeVariables.m_pathwayMiddleEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMaxEnergyMean", electronTreeVariables.m_pathwayMaxEnergyMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMinEnergySigma", electronTreeVariables.m_pathwayMinEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMiddleEnergySigma", electronTreeVariables.m_pathwayMiddleEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayMaxEnergySigma", electronTreeVariables.m_pathwayMaxEnergySigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartLength", electronTreeVariables.m_postShowerStartLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNHits", electronTreeVariables.m_postShowerStartNHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNHitsW", electronTreeVariables.m_postShowerStartNHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinNPostShowerStartHits", electronTreeVariables.m_minNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleNPostShowerStartHits", electronTreeVariables.m_middleNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxNPostShowerStartHits", electronTreeVariables.m_maxNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartScatterAngle", electronTreeVariables.m_postShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartScatterAngleW", electronTreeVariables.m_postShowerStartScatterAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartScatterAngle", electronTreeVariables.m_minPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartScatterAngle", electronTreeVariables.m_middlePostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartScatterAngle", electronTreeVariables.m_maxPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartOpeningAngleW", electronTreeVariables.m_postShowerStartOpeningAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartOpeningAngle", electronTreeVariables.m_minPostShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartOpeningAngle", electronTreeVariables.m_middlePostShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartOpeningAngle", electronTreeVariables.m_maxPostShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartOpeningAngleAsymmetryW", electronTreeVariables.m_postShowerStartOpeningAngleAsymmetryW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartOpeningAngleAsymmetry", electronTreeVariables.m_minPostShowerStartOpeningAngleAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartOpeningAngleAsymmetry", electronTreeVariables.m_middlePostShowerStartOpeningAngleAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartOpeningAngleAsymmetry", electronTreeVariables.m_maxPostShowerStartOpeningAngleAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNuVertexHitAsymmetryW", electronTreeVariables.m_postShowerStartNuVertexHitAsymmetryW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartNuVertexHitAsymmetry", electronTreeVariables.m_minPostShowerStartNuVertexHitAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartNuVertexHitAsymmetry", electronTreeVariables.m_middlePostShowerStartNuVertexHitAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexHitAsymmetry", electronTreeVariables.m_maxPostShowerStartNuVertexHitAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNuVertexEnergyAsymmetryW", electronTreeVariables.m_postShowerStartNuVertexEnergyAsymmetryW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartNuVertexEnergyAsymmetry", electronTreeVariables.m_minPostShowerStartNuVertexEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartNuVertexEnergyAsymmetry", electronTreeVariables.m_middlePostShowerStartNuVertexEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexEnergyAsymmetry", electronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartShowerStartHitAsymmetryW", electronTreeVariables.m_postShowerStartShowerStartHitAsymmetryW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartShowerStartHitAsymmetry", electronTreeVariables.m_minPostShowerStartShowerStartHitAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartShowerStartHitAsymmetry", electronTreeVariables.m_middlePostShowerStartShowerStartHitAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartShowerStartHitAsymmetry", electronTreeVariables.m_maxPostShowerStartShowerStartHitAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartShowerStartEnergyAsymmetryW", electronTreeVariables.m_postShowerStartShowerStartEnergyAsymmetryW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartShowerStartEnergyAsymmetry", electronTreeVariables.m_minPostShowerStartShowerStartEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartShowerStartEnergyAsymmetry", electronTreeVariables.m_middlePostShowerStartShowerStartEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartShowerStartEnergyAsymmetry", electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNuVertexMeanRadialDistanceW", electronTreeVariables.m_postShowerStartNuVertexMeanRadialDistanceW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartNuVertexMeanRadialDistance", electronTreeVariables.m_minPostShowerStartNuVertexMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartNuVertexMeanRadialDistance", electronTreeVariables.m_middlePostShowerStartNuVertexMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexMeanRadialDistance", electronTreeVariables.m_maxPostShowerStartNuVertexMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNuVertexEnergyWeightedMeanRadialDistanceW", electronTreeVariables.m_postShowerStartNuVertexEnergyWeightedMeanRadialDistanceW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", electronTreeVariables.m_minPostShowerStartNuVertexEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance", electronTreeVariables.m_middlePostShowerStartNuVertexEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", electronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartShowerStartMeanRadialDistanceW", electronTreeVariables.m_postShowerStartShowerStartMeanRadialDistanceW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartShowerStartMeanRadialDistance", electronTreeVariables.m_minPostShowerStartShowerStartMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartShowerStartMeanRadialDistance", electronTreeVariables.m_middlePostShowerStartShowerStartMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartShowerStartMeanRadialDistance", electronTreeVariables.m_maxPostShowerStartShowerStartMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartShowerStartEnergyWeightedMeanRadialDistanceW", electronTreeVariables.m_postShowerStartShowerStartEnergyWeightedMeanRadialDistanceW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartShowerStartEnergyWeightedMeanRadialDistance", electronTreeVariables.m_minPostShowerStartShowerStartEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance", electronTreeVariables.m_middlePostShowerStartShowerStartEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance", electronTreeVariables.m_maxPostShowerStartShowerStartEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartNuVertexMoliereRadiusW", electronTreeVariables.m_postShowerStartNuVertexMoliereRadiusW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartNuVertexMoliereRadius", electronTreeVariables.m_minPostShowerStartNuVertexMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartNuVertexMoliereRadius", electronTreeVariables.m_middlePostShowerStartNuVertexMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexMoliereRadius", electronTreeVariables.m_maxPostShowerStartNuVertexMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartShowerStartMoliereRadiusW", electronTreeVariables.m_postShowerStartShowerStartMoliereRadiusW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartShowerStartMoliereRadius", electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddlePostShowerStartShowerStartMoliereRadius", electronTreeVariables.m_middlePostShowerStartShowerStartMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartShowerStartMoliereRadius", electronTreeVariables.m_maxPostShowerStartShowerStartMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PositiveOpeningAngleW", electronTreeVariables.m_positiveOpeningAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "NegativeOpeningAngleW", electronTreeVariables.m_negativeOpeningAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxOpeningAngleW", electronTreeVariables.m_maxOpeningAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "ShowerApexLW", electronTreeVariables.m_showerApexLW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinShowerApexL", electronTreeVariables.m_minShowerApexL));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleShowerApexL", electronTreeVariables.m_middleShowerApexL));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxShowerApexL", electronTreeVariables.m_maxShowerApexL));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "ShowerApexTW", electronTreeVariables.m_showerApexTW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinShowerApexT", electronTreeVariables.m_minShowerApexT));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleShowerApexT", electronTreeVariables.m_middleShowerApexT));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxShowerApexT", electronTreeVariables.m_maxShowerApexT));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "FoundHitRatioW", electronTreeVariables.m_foundHitRatioW))
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinFoundHitRatio", electronTreeVariables.m_minFoundHitRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleFoundHitRatio", electronTreeVariables.m_middleFoundHitRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxFoundHitRatio", electronTreeVariables.m_maxFoundHitRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "FitShowerStartLW", electronTreeVariables.m_fitShowerStartLW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "FitShowerStartTW", electronTreeVariables.m_fitShowerStartTW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartMinHalfOpeningAngle", electronTreeVariables.m_postShowerStartMinHalfOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartMaxHalfOpeningAngle", electronTreeVariables.m_postShowerStartMaxHalfOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartOpeningAngle", electronTreeVariables.m_postShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartOpeningAngleAsymmetry", electronTreeVariables.m_postShowerStartOpeningAngleAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartMeanTransverseAngle", electronTreeVariables.m_postShowerStartMeanTransverseAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartMeanLWeightedTransverseAngle", electronTreeVariables.m_postShowerStartMeanLWeightedTransverseAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartMeanRadialDistance", electronTreeVariables.m_postShowerStartMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartRadialDistanceSigma", electronTreeVariables.m_postShowerStartRadialDistanceSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartEnergyWeightedMeanRadialDistance", electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartEstimatedMoliereRadius", electronTreeVariables.m_postShowerStartEstimatedMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartLWeightedMeanRadialDistance", electronTreeVariables.m_postShowerStartLWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartLWeightedRadialDistanceSigma", electronTreeVariables.m_postShowerStartLWeightedRadialDistanceSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartInitialGapSize", electronTreeVariables.m_postShowerStartInitialGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PostShowerStartMaxGapSize", electronTreeVariables.m_postShowerStartMaxGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "InitialRegionDistanceToNuVertex", electronTreeVariables.m_initialRegionDistanceToNuVertex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "InitialRegionDistanceInGaps", electronTreeVariables.m_initialRegionDistanceInGaps));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "InitialRegionMaxGapSize", electronTreeVariables.m_initialRegionMaxGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "InitialGapSizeW", electronTreeVariables.m_initialGapSizeW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinInitialGapSize", electronTreeVariables.m_minInitialGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleInitialGapSize", electronTreeVariables.m_middleInitialGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxInitialGapSize", electronTreeVariables.m_maxInitialGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "LargestGapSizeW", electronTreeVariables.m_largestGapSizeW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinLargestGapSize", electronTreeVariables.m_minLargestGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleLargestGapSize", electronTreeVariables.m_middleLargestGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxLargestGapSize", electronTreeVariables.m_maxLargestGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "LargestProjectedGapSizeW", electronTreeVariables.m_largestProjectedGapSizeW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinLargestProjectedGapSize", electronTreeVariables.m_minLargestProjectedGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleLargestProjectedGapSize", electronTreeVariables.m_middleLargestProjectedGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxLargestProjectedGapSize", electronTreeVariables.m_maxLargestProjectedGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "HitLineDensityW", electronTreeVariables.m_hitLineDensityW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinHitLineDensity", electronTreeVariables.m_minHitLineDensity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MiddleHitLineDensity", electronTreeVariables.m_middleHitLineDensity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxHitLineDensity", electronTreeVariables.m_maxHitLineDensity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "NViewsWithAmbiguousHits", electronTreeVariables.m_nViewsWithAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "NAmbiguousHits2D", electronTreeVariables.m_nAmbiguousHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinNAmbiguousHits", electronTreeVariables.m_minNAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxNAmbiguousHits", electronTreeVariables.m_maxNAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitUnaccountedEnergyU", electronTreeVariables.m_ambiguousHitUnaccountedEnergyU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitUnaccountedEnergyV", electronTreeVariables.m_ambiguousHitUnaccountedEnergyV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitUnaccountedEnergyW", electronTreeVariables.m_ambiguousHitUnaccountedEnergyW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitMinUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMinUnaccountedEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitMaxUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitShowerEnergyRatioU", electronTreeVariables.m_ambiguousHitShowerEnergyRatioU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitShowerEnergyRatioV", electronTreeVariables.m_ambiguousHitShowerEnergyRatioV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitShowerEnergyRatioW", electronTreeVariables.m_ambiguousHitShowerEnergyRatioW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitMinShowerEnergyRatio", electronTreeVariables.m_ambiguousHitMinShowerEnergyRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitMaxShowerEnergyRatio", electronTreeVariables.m_ambiguousHitMaxShowerEnergyRatio));

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

bool ShowerStartRefinementAlgorithm::TMVAIsElectron(LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables, const ParticleFlowObject *const pShowerPfo, const bool alterMetadata)
{
    m_TMVAElectronTreeVariables.m_pathwayLengthMin = electronTreeVariables.m_pathwayLengthMin;
    m_TMVAElectronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D = electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D;
    m_TMVAElectronTreeVariables.m_maxNPostShowerStartHits = electronTreeVariables.m_maxNPostShowerStartHits;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartScatterAngle = electronTreeVariables.m_maxPostShowerStartScatterAngle;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry = electronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry = electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = electronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    m_TMVAElectronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius = electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartOpeningAngle = electronTreeVariables.m_maxPostShowerStartOpeningAngle;
    m_TMVAElectronTreeVariables.m_maxFoundHitRatio = electronTreeVariables.m_maxFoundHitRatio;
    m_TMVAElectronTreeVariables.m_maxInitialGapSize = electronTreeVariables.m_maxInitialGapSize;
    m_TMVAElectronTreeVariables.m_minLargestProjectedGapSize = electronTreeVariables.m_minLargestProjectedGapSize;
    m_TMVAElectronTreeVariables.m_nViewsWithAmbiguousHits = electronTreeVariables.m_nViewsWithAmbiguousHits;
    m_TMVAElectronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy = electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy;

    float bdtScore(m_TMVAReader.EvaluateMVA("BDTG"));

    if (alterMetadata)
    {
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["ElectronConnectionPathwayScore"] = bdtScore;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pShowerPfo, metadata));
    }

    std::cout << "bdtScore: " << bdtScore << std::endl;

    return bdtScore > m_electronTMVACut;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementAlgorithm::TMVAIsGamma(LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables, const ParticleFlowObject *const pShowerPfo)
{
    m_TMVAElectronTreeVariables.m_pathwayLengthMin = electronTreeVariables.m_pathwayLengthMin;
    m_TMVAElectronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D = electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D;
    m_TMVAElectronTreeVariables.m_maxNPostShowerStartHits = electronTreeVariables.m_maxNPostShowerStartHits;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartScatterAngle = electronTreeVariables.m_maxPostShowerStartScatterAngle;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry = electronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry = electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = electronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    m_TMVAElectronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius = electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius;
    m_TMVAElectronTreeVariables.m_maxPostShowerStartOpeningAngle = electronTreeVariables.m_maxPostShowerStartOpeningAngle;
    m_TMVAElectronTreeVariables.m_maxFoundHitRatio = electronTreeVariables.m_maxFoundHitRatio;
    m_TMVAElectronTreeVariables.m_maxInitialGapSize = electronTreeVariables.m_maxInitialGapSize;
    m_TMVAElectronTreeVariables.m_minLargestProjectedGapSize = electronTreeVariables.m_minLargestProjectedGapSize;
    m_TMVAElectronTreeVariables.m_nViewsWithAmbiguousHits = electronTreeVariables.m_nViewsWithAmbiguousHits;
    m_TMVAElectronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy = electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy;

    float bdtScore(m_TMVAReader.EvaluateMVA("BDTG"));

    return bdtScore < m_gammaTMVACut;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// our signal here is actually gammas that have made a mistake by getting back to the nu vertex. i.e. their reco vertex is closer to the nu vertex than the truth says

bool ShowerStartRefinementAlgorithm::IsGamma(const ParticleFlowObject *const pPfo, const CartesianVector &nuVertexPosition) const
{
    const MCParticle *pMainMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    int pdg(std::abs(pMainMCParticle->GetParticleId()));

    if (pdg == 11)
        return false;

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

    CaloHitList caloHitList3D;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList3D);

    if (caloHitList3D.empty())
        return false;

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

    const CartesianVector &mcVertex(pMainMCParticle->GetVertex());
    const float mcSeparation((mcVertex - nuVertexPosition).GetMagnitude());
    const float recoSeparation((showerVertex - nuVertexPosition).GetMagnitude());
    const float difference(recoSeparation - mcSeparation);

    //////////////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertex, "mcVertex", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerVertex, "showerVertex", GREEN, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    //////////////////////////////////////

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

void ShowerStartRefinementAlgorithm::SetElectronTreeMetadata(const ParticleFlowObject *const pShowerPfo, 
    LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    object_creation::ParticleFlowObject::Metadata metadata;

    metadata.m_propertiesToAdd["PathwayLengthMin"] = electronTreeVariables.m_pathwayLengthMin;
    metadata.m_propertiesToAdd["MaxShowerStartPathwayScatteringAngle2D"] = electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D;
    metadata.m_propertiesToAdd["MaxNPostShowerStartHits"] = electronTreeVariables.m_maxNPostShowerStartHits;
    metadata.m_propertiesToAdd["MaxPostShowerStartScatterAngle"] = electronTreeVariables.m_maxPostShowerStartScatterAngle;
    metadata.m_propertiesToAdd["MaxPostShowerStartNuVertexEnergyAsymmetry"] = electronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry;
    metadata.m_propertiesToAdd["MaxPostShowerStartShowerStartEnergyAsymmetry"] = electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry;
    metadata.m_propertiesToAdd["MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance"] = electronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance;
    metadata.m_propertiesToAdd["MinPostShowerStartShowerStartMoliereRadius"] = electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius;
    metadata.m_propertiesToAdd["MaxPostShowerStartOpeningAngle"] = electronTreeVariables.m_maxPostShowerStartOpeningAngle;
    metadata.m_propertiesToAdd["MaxFoundHitRatio"] = electronTreeVariables.m_maxFoundHitRatio;
    metadata.m_propertiesToAdd["MaxInitialGapSize"] = electronTreeVariables.m_maxInitialGapSize;
    metadata.m_propertiesToAdd["MinLargestProjectedGapSize"] = electronTreeVariables.m_minLargestProjectedGapSize;
    metadata.m_propertiesToAdd["NViewsWithAmbiguousHits"] = electronTreeVariables.m_nViewsWithAmbiguousHits;
    metadata.m_propertiesToAdd["AmbiguousHitMaxUnaccountedEnergy"] = electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pShowerPfo, metadata));
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
        XmlHelper::ReadValue(xmlHandle, "CreateTrainingTrees", m_createTrainingTrees));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HybridMode", m_hybridMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ThresholdSignalGammaDisplacement", m_thresholdSignalGammaDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ElectronTMVACut", m_electronTMVACut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GammaTMVACut", m_gammaTMVACut));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ConnectionBDTWeightsPath", m_connectionBDTWeightsPath));

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
