/**
 *  @file   larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the svm pfo characterisation algorithm class.
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
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.h"
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
    SvmPfoCharacterisationAlgorithm::RecoParameters::RecoParameters() :
    m_minPrimaryGoodHits(15),
    m_minHitsForGoodView(5),
    m_minPrimaryGoodViews(2),
    m_foldToPrimaries(false),
    m_minHitSharingFraction(0.9f)
{
}
SvmPfoCharacterisationAlgorithm::SvmPfoCharacterisationAlgorithm() :
    m_trainingSetMode(false),
    m_enableProbability(true),
    m_useThreeDInformation(true),
    m_minProbabilityCut(0.5f),
    m_minCaloHitsCut(5),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
	m_writeToTree(false)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmPfoCharacterisationAlgorithm::Run()
{
   return PfoCharacterisationBaseAlgorithm::Run();
}


//------------------------------------------------------------------------------------------------------------------------------------------

SvmPfoCharacterisationAlgorithm::~SvmPfoCharacterisationAlgorithm()
{
  if (m_writeToTree)
    {
      //std::string m_treeName = "fufufu";
      //std::string m_fileName = "fufufu.root";
      PandoraMonitoringApi::SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool SvmPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
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
        return LArMvaHelper::Classify(m_adaBoostDecisionTree, featureVector);
    }
    else
    {
        return (LArMvaHelper::CalculateProbability(m_adaBoostDecisionTree, featureVector) > m_minProbabilityCut);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
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
	//std::cout << "testing testing testings" << std::endl;
    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    //charge related features are only calculated using hits in W view
    // This won't work unless use 3D info is set to true - dev purposes only
    const PfoCharacterisationFeatureTool::FeatureToolVector &chosenFeatureToolVector(wClusterList.empty() ? m_featureToolVectorNoChargeInfo : m_featureToolVectorThreeD);

    //std::string m_treeName = "fufufu"; // TODO This should be a member variable

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
    //std::cout << "nHits: " << nHits << std::endl;
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits", nHits);
	
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Mapping target MCParticles -> truth associated Hits
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    FillMCToRecoHitsMap(pMCParticleList, pCaloHitList, targetMCParticleToHitsMap);
    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetUnfoldedPfoToReconstructable2DHitsMap(myPfoList, targetMCParticleToHitsMap, pfoToHitsMap);
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {targetMCParticleToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
	const CaloHitList &allHitsInPfo(pfoToHitsMap.at(pPfo));
	const int nHitsInPfoTotal(allHitsInPfo.size());
	int nHitsInBestMCParticleTotal(-1), bestMCParticlePdgCode(0), bestMCParticleIsTrack(-1);
	int nHitsSharedWithBestMCParticleTotal(-1);
    CartesianVector threeDVertexPosition(0.f, 0.f, 0.f); // Mousam Vertex
    float mcEnergy = 0.f;
    float hitsShower = 0;
	float hitsTrack = 0;
	const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCHitSharingMap.at(pPfo));
	for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
	{
	    const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);
		//std::cout << "MCParticle: " << pAssociatedMCParticle->GetParticleId() << std::endl;
	    const CaloHitList &allMCHits(targetMCParticleToHitsMap.at(pAssociatedMCParticle));
		//std::cout << "allMCHits: " << allMCHits.size() << std::endl;	  
	    const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);

		if ((abs(pAssociatedMCParticle->GetParticleId()) == 11) || (pAssociatedMCParticle->GetParticleId()) == 22)
		{
			hitsShower = hitsShower + associatedMCHits.size();
		}
		else
		{
			hitsTrack = hitsTrack + associatedMCHits.size();
		}
		//std::cout << "associatedMCHits: " << associatedMCHits.size() << std::endl;
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
	//std::cout << "hitsShower: " << hitsShower << ", hitsTrack: " << hitsTrack << std::endl;
	float trackShowerHitsRatio; 
	trackShowerHitsRatio = hitsTrack/(hitsTrack + hitsShower);
	//std::cout << "trackShowerHitsRatio: " << trackShowerHitsRatio << std::endl;
	int newTrueTrackInt = (trackShowerHitsRatio >= 0.5 ? 1 : 0);
	//std::cout << "newTrueTrackInt: " << newTrueTrackInt << std::endl;
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "newTrueTrackInt", newTrueTrackInt);
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
/*
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
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionType);*/
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
	
	float trueKineticEnergy (0.f);

	if (abs(pdgCode) == 11)
		trueKineticEnergy = (mcEnergy - 0.000511);
	else if (abs(pdgCode) == 13)
		trueKineticEnergy = (mcEnergy - 0.106);
	else if (abs(pdgCode) == 2212)
		trueKineticEnergy = (mcEnergy - 0.938);
	else if (abs(pdgCode) == 211)
		trueKineticEnergy = (mcEnergy - 0.140);
	else if (abs(pdgCode) == 3312)
		trueKineticEnergy = (mcEnergy - 1.322);
	else if (abs(pdgCode) == 3222)
		trueKineticEnergy = (mcEnergy - 1.189);
	else if (abs(pdgCode) == 2112)
		trueKineticEnergy = (mcEnergy - 0.940);
	else if (abs(pdgCode) == 22)
		trueKineticEnergy = (mcEnergy);
	else if (abs(pdgCode) == 321)
		trueKineticEnergy = (mcEnergy - 0.493);
	else if (abs(pdgCode) == 311)
		trueKineticEnergy = (mcEnergy - 0.498);
	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueKineticEnergy", trueKineticEnergy);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Completeness", completeness);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Purity", purity);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pdgCode", pdgCode);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsSharedWithBestMCParticleTotal", nHitsSharedWithBestMCParticleTotal);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleTotal", nHitsInBestMCParticleTotal);
    PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoTotal", nHitsInPfoTotal);

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
    delete reader;
    // TODO Use BDT output to return from this method, so as to evaluate final output performance
    // TODO Tune this value
	//return (bdtValue > 0.0);

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

    //if no failures, proceed with svm classification
    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify((wClusterList.empty() ? m_adaBoostDecisionTreeNoChargeInfo : m_adaBoostDecisionTree), featureVector);
    }
    else
    {
        const double score(LArMvaHelper::CalculateProbability((wClusterList.empty() ? m_adaBoostDecisionTreeNoChargeInfo : m_adaBoostDecisionTree), featureVector));
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["TrackScore"] = score;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
		/*if (mcEnergy >=1.0 && mcEnergy < 1.5)
		{
		std::ofstream outfile ("sklearnBdtOutputAllPfos_1000_1500MeV.txt", std::ios::app);
		outfile << TrueTrackInt << " " << int(bool(m_minProbabilityCut <= score)) << std::endl;
		outfile.close();
		}*/
		PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "sklearnScore", score);
		int sklearnTrackInt = ( m_minProbabilityCut <= score ? 1 : 0);
		PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "sklearnTrackInt", sklearnTrackInt);
		//std::cout << "sklearnTrackInt: " << sklearnTrackInt << std::endl;
		return (m_minProbabilityCut <= score);
    }


    // End variable writing


    //return (bdtValue > 0.0);
    PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());
}

  void SvmPfoCharacterisationAlgorithm::SelectParticlesByHitCount(const MCParticleVector &candidateTargets, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap) const
{
    // Apply restrictions on the number of good hits associated with the MCParticles
    for (const MCParticle * const pMCTarget : candidateTargets)
    {
        LArMCParticleHelper::MCContributionMap::const_iterator trueHitsIter = mcToTrueHitListMap.find(pMCTarget);
        if (mcToTrueHitListMap.end() == trueHitsIter)
            continue;

        const CaloHitList &caloHitList(trueHitsIter->second);

        // Remove shared hits where target particle deposits below threshold energy fraction
        CaloHitList goodCaloHitList;
        this->SelectGoodCaloHits(&caloHitList, mcToTargetMCMap, goodCaloHitList);

        if (goodCaloHitList.size() < m_recoParameters.m_minPrimaryGoodHits)
            continue;

        unsigned int nGoodViews(0);
        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, goodCaloHitList) >= m_recoParameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, goodCaloHitList) >= m_recoParameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, goodCaloHitList) >= m_recoParameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (nGoodViews < m_recoParameters.m_minPrimaryGoodViews)
            continue;

        if (!selectedMCParticlesToHitsMap.insert(LArMCParticleHelper::MCContributionMap::value_type(pMCTarget, caloHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

void SvmPfoCharacterisationAlgorithm::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, CaloHitList &selectedGoodCaloHitList) const
{

    for (const CaloHit *const pCaloHit : *pSelectedCaloHitList)
    {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

	// fold back weights to target MC particle (primary or not)
	// keep for foldToPrimaries == false since will neglect hits belonging to MC particles
	// that have been removed at an earlier stage 
        MCParticleWeightMap weightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToTargetMCMap.find(pMCParticle);

            if (mcToTargetMCMap.end() != mcIter)
                weightMap[mcIter->second] += weight;
        }

        MCParticleVector mcTargetVector;
        for (const auto &mapEntry : weightMap) mcTargetVector.push_back(mapEntry.first);
        std::sort(mcTargetVector.begin(), mcTargetVector.end(), PointerLessThan<MCParticle>());

        const MCParticle *pBestTargetParticle(nullptr);
        float bestTargetWeight(0.f), targetWeightSum(0.f);

        for (const MCParticle *const pTargetMCParticle : mcTargetVector)
        {
            const float targetWeight(weightMap.at(pTargetMCParticle));
            targetWeightSum += targetWeight;

            if (targetWeight > bestTargetWeight)
            {
                bestTargetWeight = targetWeight;
                pBestTargetParticle = pTargetMCParticle;
            }
        }

        if (!pBestTargetParticle || (targetWeightSum < std::numeric_limits<float>::epsilon()) || ((bestTargetWeight / targetWeightSum) < m_recoParameters.m_minHitSharingFraction))
            continue;

        selectedGoodCaloHitList.push_back(pCaloHit);
    }
}

  void SvmPfoCharacterisationAlgorithm::GetMCToSelfMap(const MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToSelfMap) const
{
  for(const MCParticle *const pMCParticle : *pMCParticleList)
  {
    mcToSelfMap[pMCParticle] = pMCParticle;
  }    

}
  void SvmPfoCharacterisationAlgorithm::FillMCToRecoHitsMap(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) const
{

    // Obtain map: [MC particle -> self] (to prevent folding to primary MC particle)
    // or [MC particle -> primary mc particle] (to fold to primary MC particle)
    LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
    m_recoParameters.m_foldToPrimaries ? LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToTargetMCMap) : GetMCToSelfMap(pMCParticleList, mcToTargetMCMap);

    //REMOVED NEUTRON AND PHOTON CONSIDERATION

    // Obtain maps: [hits -> MC particle (either primary or a downstream MC particle)], [MC particle -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap trueHitToTargetMCMap;
    LArMCParticleHelper::MCContributionMap targetMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToTargetMCMap, trueHitToTargetMCMap, targetMCToTrueHitListMap);

    // Obtain vector: all or primary mc particles
    MCParticleVector targetMCVector;
    if(m_recoParameters.m_foldToPrimaries){ 
      LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, targetMCVector);
    } else {
      std::copy(pMCParticleList->begin(), pMCParticleList->end(), std::back_inserter(targetMCVector));
    }

    //REMOVED WHETHER PARTICLE MATCHES SOME CRITERIA (e.g whether downstream of neutrino) - not needed for created neutrino events

    // Remove hits that do not meet minimum hit count and share criteria
    this->SelectParticlesByHitCount(targetMCVector, targetMCToTrueHitListMap, mcToTargetMCMap, mcToRecoHitsMap);
    //SelectParticlesByHitCount(targetMCVector, targetMCToTrueHitListMap, mcToTargetMCMap, mcToRecoHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodHits", m_recoParameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForGoodView", m_recoParameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodViews", m_recoParameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_recoParameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_recoParameters.m_foldToPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseThreeDInformation", m_useThreeDInformation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmFileName", m_svmFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmName", m_svmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree)); // added by Mousam

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputTree", m_treeName)); // added by Mousam
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputFile", m_fileName)); // added by Mousam

    if (m_useThreeDInformation)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmFileNameNoChargeInfo", m_svmFileNameNoChargeInfo));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmNameNoChargeInfo", m_svmNameNoChargeInfo));
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
        if (m_svmFileName.empty() || m_svmName.empty())
        {
            std::cout << "SvmPfoCharacterisationAlgorithm: SvmFileName and SvmName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullSvmFileName(LArFileHelper::FindFileInPath(m_svmFileName, m_filePathEnvironmentVariable));
        m_adaBoostDecisionTree.Initialize(fullSvmFileName, m_svmName);

        if (m_useThreeDInformation)
        {
            if (m_svmFileNameNoChargeInfo.empty() || m_svmNameNoChargeInfo.empty())
            {
                std::cout << "SvmPfoCharacterisationAlgorithm: SvmFileName and SvmName must be set if in classification mode for no charge info in 3D mode " << std::endl;
                return STATUS_CODE_INVALID_PARAMETER;
            }
            const std::string fullSvmFileNameNoChargeInfo(LArFileHelper::FindFileInPath(m_svmFileNameNoChargeInfo, m_filePathEnvironmentVariable));
            m_adaBoostDecisionTreeNoChargeInfo.Initialize(fullSvmFileNameNoChargeInfo, m_svmNameNoChargeInfo);
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

} // namespace lar_content
