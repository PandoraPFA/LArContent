/**
 *  @file   larpandoracontent/LArMonitoring/ParticleEfficiencyAlgorithm
 *
 *  @brief  Implementation of the performance assessment algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArMonitoring/ParticleEfficiencyAlgorithm.h"

using namespace pandora;

namespace lar_content 
	{
  	ParticleEfficiencyAlgorithm::RecoParameters::RecoParameters() :
    m_minPrimaryGoodHits(15),
    m_minHitsForGoodView(5),
    m_minPrimaryGoodViews(2),
    m_foldToPrimaries(false),
    m_minHitSharingFraction(0.9f)
  	{
  	}
//------------------------------------------------------------------------------------------------------------------------------------------

  ParticleEfficiencyAlgorithm::ParticleEfficiencyAlgorithm() :
    m_caloHitListName(), 
    m_pfoListName(), 
    m_recoParameters(), 
    m_writeToTree(false), 
    m_treeName(),
    m_fileName(), 
    m_eventNumber(0)
  {
  }


  ParticleEfficiencyAlgorithm::~ParticleEfficiencyAlgorithm()
	{
    	if(m_writeToTree) 
		{
			try 
			{
        		PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
      		}
      		catch (const StatusCodeException &) 
			{
        		std::cout << "UNABLE TO WRITE TREE" << std::endl; 
      		}
    	}
  	}

  StatusCode ParticleEfficiencyAlgorithm::Run() {

    ++m_eventNumber;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
   
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Get MC Particle to reconstructable hits map
    LArMCParticleHelper::MCContributionMap mcToRecoHitsMap;
    FillMCToRecoHitsMap(pMCParticleList, pCaloHitList, mcToRecoHitsMap);

    MCParticleVector orderedTargetMCParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcToRecoHitsMap}, orderedTargetMCParticleVector);

    // If no pfos are created, still need to fill MC tree with this info
    const PfoList *pPfoList = nullptr;
    StatusCode pfoListStatus(PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    // Handle case where no pfos are created
    if(pfoListStatus == STATUS_CODE_NOT_INITIALIZED){
      AddNoPfoEntryToTree(orderedTargetMCParticleVector, mcToRecoHitsMap);
      return STATUS_CODE_NOT_INITIALIZED;
    } 
 
    // ensure that other StatusCodes are still handled
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, pfoListStatus);

    // Get pfo to reconstructable hits map
    LArMCParticleHelper::PfoContributionMap pfoToRecoHitsMap;

    if (m_recoParameters.m_foldToPrimaries) {
      // Get list of pfos to be matched with the target MC particles
      PfoList finalStatePfos;
      for (const ParticleFlowObject *const pPfo : *pPfoList) {
	if (LArPfoHelper::IsFinalState(pPfo))
	  finalStatePfos.push_back(pPfo);
      }
      LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, mcToRecoHitsMap, pfoToRecoHitsMap);
    } else {
      LArMCParticleHelper::GetUnfoldedPfoToReconstructable2DHitsMap(*pPfoList, mcToRecoHitsMap, pfoToRecoHitsMap);
    }

    PfoVector orderedPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToRecoHitsMap, orderedPfoVector);

    // Find hits that they share
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToRecoHitsMap, {mcToRecoHitsMap}, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Calculate purity and completeness for MC->Pfo matches
    LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap mcParticleToPfoCompletenessMap;
    LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap mcParticleToPfoPurityMap; 
    LArMCParticleHelper::GetMCToPfoCompletenessPurityMaps(mcToRecoHitsMap, pfoToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);


    if(m_writeToTree) {
      AddMatchesEntryToTree(orderedTargetMCParticleVector, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
    }

 
    return STATUS_CODE_SUCCESS;

  }


//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::FillMCToRecoHitsMap(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap)
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

  void ParticleEfficiencyAlgorithm::SelectParticlesByHitCount(const MCParticleVector &candidateTargets, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap)
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

void ParticleEfficiencyAlgorithm::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, CaloHitList &selectedGoodCaloHitList)
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
   
//------------------------------------------------------------------------------------------------------------------------------------------


  void ParticleEfficiencyAlgorithm::GetMCToSelfMap(const MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToSelfMap)
{
  for(const MCParticle *const pMCParticle : *pMCParticleList)
  {
    mcToSelfMap[pMCParticle] = pMCParticle;
  }    

}


//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::AddMatchesEntryToTree(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap) {

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {

      //Fill reconstructable MC particle information
      CartesianVector zUnitVector(0,0,1);
      float angleFromZ = pMCParticle->GetMomentum().GetOpeningAngle(zUnitVector);

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "EventNumber", m_eventNumber));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "MCParticleID", pMCParticle->GetParticleId()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Energy", pMCParticle->GetEnergy()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Momentum", pMCParticle->GetMomentum().GetMagnitude()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "AngleFromZ", angleFromZ));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "TotHits", static_cast<int>(mcToRecoHitsMap.at(pMCParticle).size())));

      std::vector<int> sharedHitsVector;
      std::vector<double> completenessVector;
      std::vector<double> purityVector;

      for(auto pfoSharedHitPair : mcParticleToPfoHitSharingMap.at(pMCParticle)) {

	sharedHitsVector.push_back(pfoSharedHitPair.second.size());

	for(auto pfoCompletenessPair : mcParticleToPfoCompletenessMap.at(pMCParticle)) {
	  if(pfoSharedHitPair.first == pfoCompletenessPair.first)
	    completenessVector.push_back(pfoCompletenessPair.second);
	}

	for(auto pfoPurityPair : mcParticleToPfoPurityMap.at(pMCParticle)) {
	  if(pfoSharedHitPair.first == pfoPurityPair.first)
	    purityVector.push_back(pfoPurityPair.second);
	}
      }

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "SharedHitsVector", &sharedHitsVector));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "CompletenessVector", &completenessVector));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PurityVector", &purityVector));

      PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::AddNoPfoEntryToTree(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {
    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {

      //Fill reconstructable MC particle information
      CartesianVector zUnitVector(0,0,1);
      float angleFromZ = pMCParticle->GetMomentum().GetOpeningAngle(zUnitVector);

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "EventNumber", m_eventNumber));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "MCParticleID", pMCParticle->GetParticleId()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Energy", pMCParticle->GetEnergy()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Momentum", pMCParticle->GetMomentum().GetMagnitude()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "AngleFromZ", angleFromZ));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "TotHits", static_cast<int>(mcToRecoHitsMap.at(pMCParticle).size())));

      std::vector<int> sharedHitsVector;
      std::vector<double> completenessVector;
      std::vector<double> purityVector;

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "SharedHitsVector", &sharedHitsVector));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "CompletenessVector", &completenessVector));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PurityVector", &purityVector));

      PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }
  }


//------------------------------------------------------------------------------------------------------------------------------------------
  StatusCode ParticleEfficiencyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

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
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));

  return STATUS_CODE_SUCCESS;

  }


} //namespace lar_content

