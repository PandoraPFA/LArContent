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
	m_writeToTree(true)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmPfoCharacterisationAlgorithm::Run()
{
   return PfoCharacterisationBaseAlgorithm::Run();
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

    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    //charge related features are only calculated using hits in W view
    // This won't work unless use 3D info is set to true - dev purposes only
    const PfoCharacterisationFeatureTool::FeatureToolVector &chosenFeatureToolVector(wClusterList.empty() ? m_featureToolVectorNoChargeInfo : m_featureToolVectorThreeD);

	// ATTN Assume your Pfos of interest are in a PfoList called myPfoList

    // Input lists
    const PfoList myPfoList(1, pPfo);
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
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
    CartesianVector threeDVertexPosition(0.f, 0.f, 0.f); 
	float hitsShower = 0;
	float hitsTrack = 0;

	const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCHitSharingMap.at(pPfo));
	for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
	{
	    const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);
	    const CaloHitList &allMCHits(targetMCParticleToHitsMap.at(pAssociatedMCParticle));
	    const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);

		if ((abs(pAssociatedMCParticle->GetParticleId()) == E_MINUS) || (pAssociatedMCParticle->GetParticleId()) == PHOTON)
		{
			hitsShower = hitsShower + associatedMCHits.size();
		}
		else
		{
			hitsTrack = hitsTrack + associatedMCHits.size();
		}
	    if (static_cast<int>(associatedMCHits.size()) > nHitsSharedWithBestMCParticleTotal)
	    {
		 nHitsSharedWithBestMCParticleTotal = associatedMCHits.size();

		 nHitsInBestMCParticleTotal = allMCHits.size();

		 bestMCParticlePdgCode = pAssociatedMCParticle->GetParticleId();               
		 bestMCParticleIsTrack = ((PHOTON != pAssociatedMCParticle->GetParticleId()) && (E_MINUS != std::abs(pAssociatedMCParticle->GetParticleId())) ? 1 : 0);
		 threeDVertexPosition = pAssociatedMCParticle->GetVertex();
	    }
	}

	const float completeness((nHitsInBestMCParticleTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInBestMCParticleTotal) : 0.f);
    const float purity((nHitsInPfoTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInPfoTotal) : 0.f);
    float xVertexPos = threeDVertexPosition.GetX();
    float yVertexPos = threeDVertexPosition.GetY();
    float zVertexPos = threeDVertexPosition.GetZ();
//-----------------------------------------------------------------------------------------------------------
   
    const int trueTrackInt = bestMCParticleIsTrack;
 	CaloHitList threeDCaloHitList;
  	LArPfoHelper::GetCaloHits(pPfo, TPC_3D, threeDCaloHitList);
	CaloHitList checkHitListAll;

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

	mischaracterisedPfo = ((((showerProbability < 0.5) && (trueTrackInt == 0)) || ((showerProbability > 0.5) && (trueTrackInt == 1))) ? 1 : 0);

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Start variable writing
    const LArMvaHelper::MvaFeatureVector featureVector(LArMvaHelper::CalculateFeatures(chosenFeatureToolVector, this, pPfo));
    
    if (m_trainingSetMode)
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
            LArMvaHelper::ProduceTrainingExample(outputFile, isTrueTrack, featureVector);
			}
        }
        return isTrueTrack;
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
		return (m_minProbabilityCut <= score);
    }
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
