/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the shower start refinement algorithm class.
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

using namespace pandora;

namespace lar_content
{

ShowerStartRefinementAlgorithm::ShowerStartRefinementAlgorithm() : 
    m_binSize(0.005),     
    m_createTrainingTrees(false),
    m_minElectronCompleteness(0.33f),
    m_minElectronPurity(0.5f)
    //m_electronTMVACut(-0.1),
    //m_gammaTMVACut(-0.1),
    //m_TMVAReader("")
{
    /*
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
    */
  ///////
  /*
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
    //std::string weightFilePath = "TMVAClassification_BDTG_FINAL.weights.xml";
    m_TMVAReader.BookMVA("BDTG", weightFilePath);
  */
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
        }

        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ElectronBackgroundTree", "ConnectionPathwayTrees.root", "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::Run()
{
    PfoVector pfoVector;
    this->FillPfoVector(pfoVector);

    if (pfoVector.empty())
        return STATUS_CODE_SUCCESS;

    CartesianVector nuVertexPosition(0.f, 0.f, 0.f);
    if (this->GetNeutrinoVertex(nuVertexPosition) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListU", pCaloHitListU));

    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListV", pCaloHitListV));

    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitListW));

    if (m_createTrainingTrees)
        this->FillElectronHitMap();

    // Run shower refinement tools
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
// Random functions?
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
/*
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
*/
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

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pShowerPfo, metadata));
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
        XmlHelper::ReadValue(xmlHandle, "CreateTrainingTrees", m_createTrainingTrees));

    /*
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ElectronTMVACut", m_electronTMVACut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GammaTMVACut", m_gammaTMVACut));
    */

    PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

// Functions for training the BDT
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

void ShowerStartRefinementAlgorithm::FillTree(const std::string &treeName, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinLargestProjectedGapSize", electronTreeVariables.m_minLargestProjectedGapSize));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxInitialGapSize", electronTreeVariables.m_maxInitialGapSize));]
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "PathwayLengthMin", electronTreeVariables.m_pathwayLengthMin));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxShowerStartPathwayScatteringAngle2D", electronTreeVariables.m_maxShowerStartPathwayScatteringAngle2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxNPostShowerStartHits", electronTreeVariables.m_maxNPostShowerStartHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartOpeningAngle", electronTreeVariables.m_maxPostShowerStartOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartShowerStartEnergyAsymmetry", electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MinPostShowerStartShowerStartMoliereRadius", electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxFoundHitRatio", electronTreeVariables.m_maxFoundHitRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartScatterAngle", electronTreeVariables.m_maxPostShowerStartScatterAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexEnergyAsymmetry", electronTreeVariables.m_maxPostShowerStartNuVertexEnergyAsymmetry));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance", electronTreeVariables.m_maxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "NViewsWithAmbiguousHits", electronTreeVariables.m_nViewsWithAmbiguousHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "AmbiguousHitMaxUnaccountedEnergy", electronTreeVariables.m_ambiguousHitMaxUnaccountedEnergy));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName));
}

} // namespace lar_content




