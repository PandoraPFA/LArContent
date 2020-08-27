/**
 *  @file   larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include <vector>
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "Helpers/MCParticleHelper.h"
#include "larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <fstream>
#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/CutPfoCharacterisationAlgorithm.h"
using namespace pandora;

namespace lar_content
{

template<typename T>
MvaPfoCharacterisationAlgorithm<T>::MvaPfoCharacterisationAlgorithm() :
    m_trainingSetMode(false),
    m_enableProbability(true),
    m_useThreeDInformation(true),
    m_minProbabilityCut(0.5f),
    m_minCaloHitsCut(5),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_writeToTree(true),
    m_eventNumber(-1)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmPfoCharacterisationAlgorithm::Run()
{
   ++m_eventNumber;
   return PfoCharacterisationBaseAlgorithm::Run();
}


//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
MvaPfoCharacterisationAlgorithm<T>::~MvaPfoCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        std::string m_treeName = "new_BdtVertex_cp80_2hierarchy_minus_two";
        std::string m_fileName = "new_BdtVertex_cp80_2hierarchy_minus_two.root";
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool MvaPfoCharacterisationAlgorithm<T>::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    const LArMvaHelper::MvaFeatureVector featureVector(LArMvaHelper::CalculateFeatures(m_featureToolVector, this, pCluster));

    if (m_trainingSetMode)
    {
        bool isTrueTrack(false);

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
        }
        catch (const StatusCodeException &) {}
        LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, isTrueTrack, featureVector);
        return isTrueTrack;
    }

    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify(m_mva, featureVector);
    }
    else
    {
        return (LArMvaHelper::CalculateProbability(m_mva, featureVector) > m_minProbabilityCut);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool MvaPfoCharacterisationAlgorithm<T>::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    static unsigned int pfoNumber = -1;
    ++pfoNumber; 
    if (!LArPfoHelper::IsThreeD(pPfo))
    {
        if (m_enableProbability)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["TrackScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        }
        return (pPfo->GetParticleId() == MU_MINUS);
    }

    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    //charge related features are only calculated using hits in W view
    // This won't work unless use 3D info is set to true - dev purposes only
    const PfoCharacterisationFeatureTool::FeatureToolVector &chosenFeatureToolVector(wClusterList.empty() ? m_featureToolVectorNoChargeInfo : m_featureToolVectorThreeD);

    std::string m_treeName = "jhanzeb"; // TODO This should be a member variable

    // Purity, completeness
    // ATTN Assume your Pfos of interest are in a PfoList called myPfoList

    // Input lists
    const PfoList myPfoList(1, pPfo);

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

	//------------------------------------------------------------------mc neutrino energy and nHits----------------------------------------------------------------------------------------------

	MCParticleVector mcNeutrinoVector;
	LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);
	if (mcNeutrinoVector.size() > 1)
		{
		throw StatusCodeException(STATUS_CODE_FAILURE);
		}
	else
		{
		float mcNeutrinoEnergy = (mcNeutrinoVector.front())->GetEnergy();
    	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNeutrinoEnergy", mcNeutrinoEnergy); 
		}

  	CaloHitList threeDCaloHitList;
  	LArPfoHelper::GetCaloHits(pPfo, TPC_3D, threeDCaloHitList);
	int nHits = threeDCaloHitList.size();
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits", nHits);
	
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Mapping target MCParticles -> truth associated Hits
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::PrimaryParameters mcPrimaryParameters;
    mcPrimaryParameters.m_selectInputHits = false;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, mcPrimaryParameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(myPfoList, targetMCParticleToHitsMap, pfoToHitsMap);

    // Last step
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {targetMCParticleToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);

	const CaloHitList &allHitsInPfo(pfoToHitsMap.at(pPfo));
	
	const int nHitsInPfoTotal(allHitsInPfo.size());

	int nHitsInBestMCParticleTotal(-1), bestMCParticlePdgCode(0), bestMCParticleIsTrack(-1);
	int nHitsSharedWithBestMCParticleTotal(-1);
        CartesianVector threeDVertexPosition(0.f, 0.f, 0.f); // Mousam Vertex
        float mcEnergy = 0.f;
	const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCHitSharingMap.at(pPfo));

	for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
	{
	    const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);

	    const CaloHitList &allMCHits(targetMCParticleToHitsMap.at(pAssociatedMCParticle));
	  
	    const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);

	    if (static_cast<int>(associatedMCHits.size()) > nHitsSharedWithBestMCParticleTotal)
	    {
		 nHitsSharedWithBestMCParticleTotal = associatedMCHits.size();

		 nHitsInBestMCParticleTotal = allMCHits.size();

		 bestMCParticlePdgCode = pAssociatedMCParticle->GetParticleId();               
		 bestMCParticleIsTrack = ((PHOTON != pAssociatedMCParticle->GetParticleId()) && (E_MINUS != std::abs(pAssociatedMCParticle->GetParticleId())) ? 1 : 0);
		 threeDVertexPosition = pAssociatedMCParticle->GetVertex(); // Mousam Vertex
         mcEnergy = pAssociatedMCParticle->GetEnergy();
	    }
	}

    //---------------------------------------------Get Pfo Energy---------------------------------------------------   
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEnergy", mcEnergy);  

    //----------------------------------------------Get Momentum------------------------------------------------
    
    CartesianVector momentumVector(0.f, 0.f, 0.f);
    
    momentumVector = pPfo->GetMomentum();

    float momentumX = momentumVector.GetX();
    float momentumY = momentumVector.GetY();
    float momentumZ = momentumVector.GetZ();

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoMomentumX", momentumX);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoMomentumY", momentumY);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoMomentumZ", momentumZ);    
        
    //-------------------------------Interaction-Type-writing---------------------------------------------------------------------------------

    std::string interType; 

    MCParticleVector mcPrimaryVector;
    MCParticleList mcPrimaryList;
    for (const auto &mapEntry : targetMCParticleToHitsMap) mcPrimaryVector.push_back(mapEntry.first);
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
    mcPrimaryList.push_back(pMCPrimary);
    }
   
    const LArInteractionTypeHelper::InteractionType interactionTypeStr(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    
    int interactionType = static_cast<int>(interactionTypeStr);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionType);

    //------------------------------Vertex-----------------------------------------------------------------------
    float xVertexPos = threeDVertexPosition.GetX();
    float yVertexPos = threeDVertexPosition.GetY();
    float zVertexPos = threeDVertexPosition.GetZ();

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xVertexPos", xVertexPos);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "yVertexPos", yVertexPos);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zVertexPos", zVertexPos);

    //----------------------------isInFiducialVolume------------------------------------------------------------

    int isInFiducialVolume(0);

    isInFiducialVolume = (((((xVertexPos <= -20) && (xVertexPos >= -340)) || ((xVertexPos <= 340) && (xVertexPos >= 20))) && ((abs(yVertexPos) <= 584)) && (zVertexPos >= 200 || zVertexPos <= 1194)) ? 1 : 0);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isInFiducialVolume", isInFiducialVolume);
    
    //-----------------------------------------------------------------------------------------------------------

    const float completeness((nHitsInBestMCParticleTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInBestMCParticleTotal) : 0.f);
    const float purity((nHitsInPfoTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInPfoTotal) : 0.f);
    int pdgCode = bestMCParticlePdgCode;

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Completeness", completeness);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Purity", purity);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pdgCode", pdgCode);

    const int TrueTrackInt = bestMCParticleIsTrack;

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrueTrackInt", TrueTrackInt);

    // End purity, completeness

	//----------------------------veto second version----------------------------------------------------

	CaloHitList checkHitListW;
	CaloHitList checkHitListU;
	CaloHitList checkHitListV;
	CaloHitList checkHitListAll;

	LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, checkHitListW);
	LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, checkHitListU);
	LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, checkHitListV);

	checkHitListAll.splice(checkHitListAll.end(), checkHitListW);
	checkHitListAll.splice(checkHitListAll.end(), checkHitListU);
	checkHitListAll.splice(checkHitListAll.end(), checkHitListV);

	LArMCParticleHelper::MCRelationMap mcPrimaryMap;
	LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
	LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
	LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

	LArMCParticleHelper::GetMCParticleToCaloHitMatches(&checkHitListAll, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

	int showerCount(0);
	int mischaracterisedPfo(0);

	for (const CaloHit *pHit : checkHitListAll)
	{

		const MCParticle *pHitMCParticle(nullptr);

		try
		{
			pHitMCParticle = hitToMCMap.at(pHit);
		}
		catch (...) {continue;}

		if ((PHOTON == pHitMCParticle->GetParticleId()) || (E_MINUS == std::abs(pHitMCParticle->GetParticleId())))
		{
			++showerCount;
		}
	}
	
	float showerProbability = (static_cast<float>(showerCount))/(static_cast<float>(hitToMCMap.size()));

	mischaracterisedPfo = ((((showerProbability < 0.5) && (TrueTrackInt == 0)) || ((showerProbability > 0.5) && (TrueTrackInt == 1))) ? 1 : 0);

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mischaracterisedPfo", mischaracterisedPfo);
    //--------------------------------------Event and Pfo number writing-(uncomment)-------------------------------------------------------------------------

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber);

    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoNumber", static_cast<int>(pfoNumber));

    //-----------------------------------Writing 3D hit position in tree - Mousam----------------------------------------------------------------
    
    CaloHitList uCaloList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, uCaloList);
    std::vector<float> uViewXHitVector;
    std::vector<float> uViewYHitVector;
    std::vector<float> uViewZHitVector;
    
    int nUEvent(0);
    for (const CaloHit *const pCaloHit : uCaloList)
      {

	float x = pCaloHit->GetPositionVector().GetX();
    float y = pCaloHit->GetPositionVector().GetY();
	float z = pCaloHit->GetPositionVector().GetZ();
	
	uViewXHitVector.push_back(x);
	uViewYHitVector.push_back(y);        
	uViewZHitVector.push_back(z);
	++nUEvent;
      };
    
    CaloHitList vCaloList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, vCaloList);
    std::vector<float> vViewXHitVector;
    std::vector<float> vViewYHitVector;
    std::vector<float> vViewZHitVector;

    int nVEvent(0);
    for (const CaloHit *const pCaloHit : vCaloList)
      {
	float x = pCaloHit->GetPositionVector().GetX();
    float y = pCaloHit->GetPositionVector().GetY();
	float z = pCaloHit->GetPositionVector().GetZ();
	
	vViewXHitVector.push_back(x);
	vViewYHitVector.push_back(y);
	vViewZHitVector.push_back(z);
	++nVEvent;
      };
    
    CaloHitList wCaloList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wCaloList);
    std::vector<float> wViewXHitVector;
    std::vector<float> wViewYHitVector;
    std::vector<float> wViewZHitVector;

    int nWEvent(0);
    for (const CaloHit *const pCaloHit : wCaloList)
      {
	float x = pCaloHit->GetPositionVector().GetX();
    float y = pCaloHit->GetPositionVector().GetY();
	float z = pCaloHit->GetPositionVector().GetZ();
	
	wViewXHitVector.push_back(x);
	wViewYHitVector.push_back(y);
    wViewZHitVector.push_back(z);
	++nWEvent;
      };


    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "uViewXHitVector", &uViewXHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "uViewYHitVector", &uViewYHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "uViewZHitVector", &uViewZHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vViewXHitVector", &vViewXHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vViewYHitVector", &uViewYHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vViewZHitVector", &vViewZHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "wViewXHitVector", &wViewXHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "wViewYHitVector", &uViewYHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "wViewZHitVector", &wViewZHitVector);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nUEvent", nUEvent);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nVEvent", nVEvent);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nWEvent", nWEvent);

    //----------------------------------------current Pfo Characterisation------------------------------------------------------------------------
	
	float straightLineLength(-1.f);
    float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());
	bool checkVar(true);
    double m_slidingFitWindow(5);
    double m_slidingShowerFitWindow(10);
    double m_maxShowerLengthCut(500.f);
    double m_dTdLWidthRatioCut(0.08f);
    double m_vertexDistanceRatioCut(500.f);
    double m_showerWidthRatioCut(0.2f);


    ClusterList twoDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);
    int currentTrackInt = 0;
    typedef std::set<pandora::HitType> HitTypeSet;
    HitTypeSet hitTypeSet;

    unsigned int nTrackLikeViews(0);
    for (const Cluster *const pCluster : twoDClusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        if (!hitTypeSet.insert(hitType).second)
            continue;
  		try
     	{
        	const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        	straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude();

        	for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
       		{
            	dTdLMin = std::min(dTdLMin, static_cast<float>(mapEntry.second.GetGradient()));
            	dTdLMax = std::max(dTdLMax, static_cast<float>(mapEntry.second.GetGradient()));
       		}
    	}
    	catch (const StatusCodeException &)
    	{
    	}

    	if (straightLineLength < std::numeric_limits<float>::epsilon())
        	checkVar = false;

    	if (straightLineLength > m_maxShowerLengthCut)
        	checkVar = true;

    	if ((dTdLMax - dTdLMin) / straightLineLength > m_dTdLWidthRatioCut)
        	checkVar = false;

    	const float vertexDistance(CutClusterCharacterisationAlgorithm::GetVertexDistance(this, pCluster));

    	if ((vertexDistance > std::numeric_limits<float>::epsilon()) && ((vertexDistance / straightLineLength) > m_vertexDistanceRatioCut))
        	checkVar = false;

    	const float showerFitWidth(CutClusterCharacterisationAlgorithm::GetShowerFitWidth(this, pCluster, m_slidingShowerFitWindow));

    	if ((showerFitWidth < std::numeric_limits<float>::epsilon()) || ((showerFitWidth / straightLineLength) > m_showerWidthRatioCut))
        	checkVar = false;

        if (checkVar)
            ++nTrackLikeViews;

        if (nTrackLikeViews >= m_minTrackLikeViews)
            currentTrackInt = 1;
    }
	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "currentTrackInt", currentTrackInt);
    //--------------------------------------------------------------------------------------------------------------------------------------
    // Start variable writing
    const LArMvaHelper::MvaFeatureVector featureVector(LArMvaHelper::CalculateFeatures(chosenFeatureToolVector, this, pPfo));

    const LArMvaHelper::MvaFeatureVector threeDLinearFitFeatureVectorOfType(LArMvaHelper::CalculateFeaturesOfType<ThreeDLinearFitFeatureTool>(chosenFeatureToolVector, this, pPfo));
    
    float length(threeDLinearFitFeatureVectorOfType.at(0).IsInitialized() ? threeDLinearFitFeatureVectorOfType.at(0).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "length", length);

    float diff(threeDLinearFitFeatureVectorOfType.at(1).IsInitialized() ? threeDLinearFitFeatureVectorOfType.at(1).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "diff", diff);

    float gap(threeDLinearFitFeatureVectorOfType.at(2).IsInitialized() ? threeDLinearFitFeatureVectorOfType.at(2).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "gap", gap);
 	

    float rms(threeDLinearFitFeatureVectorOfType.at(3).IsInitialized() ? threeDLinearFitFeatureVectorOfType.at(3).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rms", rms);

    const LArMvaHelper::MvaFeatureVector threeDVertexDistanceFeatureVectorOfType(LArMvaHelper::CalculateFeaturesOfType<ThreeDVertexDistanceFeatureTool>(chosenFeatureToolVector, this, pPfo));
    
    float vertexDistance(threeDVertexDistanceFeatureVectorOfType.at(0).IsInitialized() ? threeDVertexDistanceFeatureVectorOfType.at(0).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance", vertexDistance);

    const LArMvaHelper::MvaFeatureVector threeDOpeningAngleFeatureVectorOfType(LArMvaHelper::CalculateFeaturesOfType<ThreeDOpeningAngleFeatureTool>(chosenFeatureToolVector, this, pPfo));
    
    float diffAngle(threeDOpeningAngleFeatureVectorOfType.at(0).IsInitialized() ? threeDOpeningAngleFeatureVectorOfType.at(0).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "diffAngle", diffAngle);

    const LArMvaHelper::MvaFeatureVector threeDPCAFeatureVectorOfType(LArMvaHelper::CalculateFeaturesOfType<ThreeDPCAFeatureTool>(chosenFeatureToolVector, this, pPfo));
    
    float pca1(threeDPCAFeatureVectorOfType.at(0).IsInitialized() ? threeDPCAFeatureVectorOfType.at(0).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pca1", pca1);

    float pca2(threeDPCAFeatureVectorOfType.at(1).IsInitialized() ? threeDPCAFeatureVectorOfType.at(1).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pca2", pca2);

    int isChargeInfoAvailable(!wClusterList.empty());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isChargeInfoAvailable", isChargeInfoAvailable);

    float charge1(-std::numeric_limits<float>::max()), charge2(-std::numeric_limits<float>::max());

    if (!wClusterList.empty())
    {
        const LArMvaHelper::MvaFeatureVector threeDChargeFeatureVectorOfType(LArMvaHelper::CalculateFeaturesOfType<ThreeDChargeFeatureTool>(chosenFeatureToolVector, this, pPfo));
        charge1 = (threeDChargeFeatureVectorOfType.at(0).IsInitialized() ? threeDChargeFeatureVectorOfType.at(0).Get() : -std::numeric_limits<float>::max());
        charge2 = (threeDChargeFeatureVectorOfType.at(1).IsInitialized() ? threeDChargeFeatureVectorOfType.at(1).Get() : -std::numeric_limits<float>::max());
    }
    
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "charge1", charge1);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "charge2", charge2);
    
    const LArMvaHelper::MvaFeatureVector pfoHierarchyFeatureVectorOfType(LArMvaHelper::CalculateFeaturesOfType<PfoHierarchyFeatureTool>(chosenFeatureToolVector, this, pPfo));
	
	float nAllDaughter(pfoHierarchyFeatureVectorOfType.at(0).IsInitialized() ? pfoHierarchyFeatureVectorOfType.at(0).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nAllDaughter", nAllDaughter);

	float nHits3DDaughterTotal(pfoHierarchyFeatureVectorOfType.at(1).IsInitialized() ? pfoHierarchyFeatureVectorOfType.at(1).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits3DDaughterTotal", nHits3DDaughterTotal);

	float daughterParentNhitsRatio(pfoHierarchyFeatureVectorOfType.at(2).IsInitialized() ? pfoHierarchyFeatureVectorOfType.at(2).Get() : -std::numeric_limits<float>::max());
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "daughterParentNhitsRatio", daughterParentNhitsRatio);

    PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());

	//std::cout << "happy" << std::endl;

    // TODO Apply trained BDT output
    // TODO FInalise variables, training, etc.
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    reader->AddVariable("length", &length);
    reader->AddVariable("diff", &diff);
    reader->AddVariable("gap", &gap);
    reader->AddVariable("rms", &rms);
    reader->AddVariable("vertexDistance", &vertexDistance);
    reader->AddVariable("diffAngle", &diffAngle);
    reader->AddVariable("pca1", &pca1);
    reader->AddVariable("pca2", &pca2);
    reader->AddVariable("charge1", &charge1);
    reader->AddVariable("charge2", &charge2);
    reader->AddVariable("nAllDaughter", &nAllDaughter);
    reader->AddVariable("nHits3DDaughterTotal", &nHits3DDaughterTotal);
    reader->AddVariable("daughterParentNhitsRatio", &daughterParentNhitsRatio);

    reader->BookMVA("BDT method", "/storage/epp2/phrwdg/Dune/newVarPandora/tmva/dataset/weights/new_BdtVertex_bdt3_type2_vertex_BDT.weights.xml");
    const float bdtValue(reader->EvaluateMVA("BDT method"));
	int bdtTrackInt = (bdtValue >= 0.0 ? 1 : 0);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bdtValue", bdtValue);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bdtTrackInt", bdtTrackInt);
    PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());
    delete reader;
    // TODO Use BDT output to return from this method, so as to evaluate final output performance
    // TODO Tune this value
	//return (bdtValue > 0.0);
    // Display interesting pfo
    /*if (length > 148 && length < 149 && abs(pdgCode) == 2212 && gap > 0.1441 && gap < 0.1442 && bdtValue > 0.071 && bdtValue < 0.072)
    {
      std::cout << "The length value is: " << length << std::endl;
      std::cout << "The vertexDistance value is: " << vertexDistance << std::endl;
      std::cout << "The conc2 value is: " << conc2 << std::endl;
      std::cout << "The conic value is: " << conic << std::endl;
      std::cout << "The densY value is: " << densY << std::endl;
      std::cout << "The densZ value is: " << densZ << std::endl;
      std::cout << "The densYRatio value is: " << densYRatio << std::endl;
      std::cout << "The densZRatio value is: " << densZRatio << std::endl;
      std::cout << "The nSegDoubleHits value is: " << nSegDoubleHits << std::endl;
      std::cout << "The nSegDoubleHitsRatio value is: " << nSegDoubleHitsRatio << std::endl;
      const PfoList myPfoList1(1, pPfo);
        PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &myPfoList1, "MyPfoList", RED, true, false);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }*/
    // End display interesting pfo

    if (m_trainingSetMode)
    {
/*      bool isTrueTrack(false);
        bool isMainMCParticleSet(false);

        try
        {
            const MCParticle *const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
            isMainMCParticleSet = (pMCParticle->GetParticleId() != 0);Completeness >= comp && Purity >= pure && (abs(xVertexPos) <= 340) && (abs(yVertexPos) <= 584) && (zVertexPos >= 200 && zVertexPos <= 1194)
        }
        catch (const StatusCodeException &) {}
*/
		{
        const bool isTrueTrack(1 == bestMCParticleIsTrack);
        const bool isMainMCParticleSet(0 != bestMCParticlePdgCode);

        if (isMainMCParticleSet)
        {
			if (completeness >= 0.8 && purity >= 0.8 && mischaracterisedPfo == 0 && (abs(xVertexPos) <= 340) && (abs(yVertexPos) <= 584) && (zVertexPos >= 200 && zVertexPos <= 1194))
			{
            std::string outputFile;
            outputFile.append(m_trainingOutputFile);
            const std::string end=((wClusterList.empty()) ? "noChargeInfo.txt" : ".txt");
            outputFile.append(end);
            LArMvaHelper::ProduceTrainingExample(outputFile, isTrueTrack, featureVector); // TODO Need this for sklearn training
			}
        }

        return isTrueTrack;
		}
    } // training mode


    for (const LArMvaHelper::MvaFeature &featureValue : featureVector)
    {
        if (!featureValue.IsInitialized())
        {
            if (m_enableProbability)
            {
                object_creation::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["TrackScore"] = -1.f;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
            return (pPfo->GetParticleId() == MU_MINUS);
        }
    }

    //if no failures, proceed with MvaPfoCharacterisationAlgorithm classification
    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify((wClusterList.empty() ? m_mvaNoChargeInfo : m_mva), featureVector);
    }
    else
    {
        const double score(LArMvaHelper::CalculateProbability((wClusterList.empty() ? m_mvaNoChargeInfo : m_mva), featureVector));
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["TrackScore"] = score;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
		//std::cout << "start!!" << std::endl;
		/*if (mcEnergy >=1.0 && mcEnergy < 1.5)
		{
		std::cout << "mcEnergy: " << mcEnergy << std::endl;
		std::ofstream outfile ("sklearnBdtOutputAllPfos_1000_1500MeV.txt", std::ios::app);
		outfile << TrueTrackInt << " " << int(bool(m_minProbabilityCut <= score)) << std::endl;
		std::cout << " hello!!!" << std::endl;
		outfile.close();
		}*/
		//std::cout << "question" << std::endl;
		PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "sklearnScore", score);
		int sklearnTrackInt = ( m_minProbabilityCut <= score ? 1 : 0);
		PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "sklearnTrackInt", sklearnTrackInt);
		//std::cout << "answer" << std::endl;
		return (m_minProbabilityCut <= score);
    }


    // End variable writing


    //return (bdtValue > 0.0);

}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode MvaPfoCharacterisationAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseThreeDInformation", m_useThreeDInformation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    // ATTN Support legacy XML configurations (note an order of precedence of XML keys exists)
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtFileName", m_mvaFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmFileName", m_mvaFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileName", m_mvaFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaName", m_mvaName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree)); // added by Mousam

    /*PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName)); // added by Mousam
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName)); // added by Mousam*/

    if (m_useThreeDInformation)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtNameNoChargeInfo", m_mvaNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmNameNoChargeInfo", m_mvaNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaNameNoChargeInfo", m_mvaNameNoChargeInfo));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableProbability",  m_enableProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinProbabilityCut", m_minProbabilityCut));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    if (m_trainingSetMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        if (m_mvaFileName.empty() || m_mvaName.empty())
        {
            std::cout << "MvaPfoCharacterisationAlgorithm: MvaFileName and MvaName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullMvaFileName(LArFileHelper::FindFileInPath(m_mvaFileName, m_filePathEnvironmentVariable));
        m_mva.Initialize(fullMvaFileName, m_mvaName);

        if (m_useThreeDInformation)
        {
            if (m_mvaFileNameNoChargeInfo.empty() || m_mvaNameNoChargeInfo.empty())
            {
                std::cout << "MvaPfoCharacterisationAlgorithm: MvaFileNameNoChargeInfo and MvaNameNoChargeInfo must be set if in classification mode for no charge info in 3D mode " << std::endl;
                return STATUS_CODE_INVALID_PARAMETER;
            }
            const std::string fullMvaFileNameNoChargeInfo(LArFileHelper::FindFileInPath(m_mvaFileNameNoChargeInfo, m_filePathEnvironmentVariable));
            m_mvaNoChargeInfo.Initialize(fullMvaFileNameNoChargeInfo, m_mvaNameNoChargeInfo);
        }
    }

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    if (m_useThreeDInformation)
    {
        AlgorithmToolVector algorithmToolVectorNoChargeInfo;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureToolsNoChargeInfo", algorithmToolVectorNoChargeInfo));
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVectorThreeD));
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVectorNoChargeInfo)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVectorNoChargeInfo));
    }
    else
    {
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));
    }

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

template class MvaPfoCharacterisationAlgorithm<AdaBoostDecisionTree>;
template class MvaPfoCharacterisationAlgorithm<SupportVectorMachine>;

} // namespace lar_content
