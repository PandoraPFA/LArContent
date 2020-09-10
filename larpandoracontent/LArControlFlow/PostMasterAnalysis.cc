//check out final cosmic/neutrino tags

#include "larpandoracontent/LArControlFlow/PostMasterAnalysis.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "PandoraMonitoringApi.h"
//#include "Api/PandoraApi.h"
//#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
//#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

using namespace pandora;

namespace lar_content
{

  PostMasterAnalysis::PostMasterAnalysis() {
  }

  PostMasterAnalysis::~PostMasterAnalysis() {
    try
      {
  	PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreec", "outputc.root", "UPDATE"));
      }
    catch (const StatusCodeException &)
      {
  	std::cout << " Unable to write tree  to file " << std::endl;
      }
  }
  //-----------------------------------------------------------------------

  StatusCode PostMasterAnalysis::Run() {

    const PfoList *pPfoList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));  //RecreatedPfos

    for (const ParticleFlowObject *const pPfo : *pPfoList)
      {
	CaloHitList totalcalohits2;
	//	std::cout << "PMA: pPfo " << pPfo << "   "  << pPfo->GetParticleId()  << std::endl;
	LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, totalcalohits2);
	LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, totalcalohits2);
	LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, totalcalohits2);
	std::cout << totalcalohits2.size() << std::endl;
	if( pPfo->GetParticleId()  ==14 ) {
	  // std::cout << "NEUTRINO!" << std::endl;
	}
      }

    int Size = 0;
    std::list<std::pair<const pandora::MCParticle*, int>> mchitslist;
    std::list<std::pair<const pandora::MCParticle*, int>> mcpairslist;
    MCParticleList mcparticlelist;
    std::list<std::pair<const pandora::MCParticle*, int>>::iterator it;
    std::list<std::pair<const pandora::MCParticle*, int>>::iterator it2;
    CaloHitList totalcalohits;
    int clearcosmiccount = 0;
    int clearneutrinocount = 0;
    int correctcr = 0;
    int incorrectcr = 0;
    int correctneutrino = 0;
    int incorrectneutrino = 0;
    int taggedneutrino = 0;
    int wasright = 0;
	//---------------------------------------------------------------------------------------------------------

    for (const Pfo *const pPPPPfo :  *pPfoList ) {
      try{
	const MCParticle *pMCParticle = LArMCParticleHelper::GetMainMCParticle(pPPPPfo);
	mcparticlelist.push_back(pMCParticle);
	LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_U, totalcalohits);
	LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_V, totalcalohits);
	LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_W, totalcalohits);

      } catch(...){
	std::cout << "lost one mc particle : " << pPPPPfo << std::endl;
      }

    }

    //get the MC particles to hits map
    LArMCParticleHelper::MCContributionMap mcMCParticlesToGoodHitsMap;
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 0;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_minPrimaryGoodViews = 0;
    parameters.m_minHitSharingFraction = 0;
    parameters.m_selectInputHits = true;
    LArMCParticleHelper::SelectReconstructableMCParticles(&mcparticlelist, &totalcalohits, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcMCParticlesToGoodHitsMap); 


    //get the pfos to hits map
    LArMCParticleHelper::PfoContributionMap PfosToGoodHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, {mcMCParticlesToGoodHitsMap}, PfosToGoodHitsMap);  //changed

    //find hits shared between Pfos and MCParticles
    LArMCParticleHelper::PfoToMCParticleHitSharingMap PfotoMCParticleMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap ParticletoPfoMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(PfosToGoodHitsMap, {mcMCParticlesToGoodHitsMap}, PfotoMCParticleMap, ParticletoPfoMap);
    
    MCParticleVector mcParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcMCParticlesToGoodHitsMap}, mcParticleVector);
    PfoVector pfoVector;
    LArMonitoringHelper::GetOrderedPfoVector({PfosToGoodHitsMap}, pfoVector);

    unsigned int nMatches = std::numeric_limits<unsigned int>::max();
    LArMonitoringHelper::PrintMatchingTable(pfoVector, mcParticleVector, ParticletoPfoMap, nMatches);

    
    for (const auto &pMCParticle : mcParticleVector)
      {
	const auto &caloHitList2 = mcMCParticlesToGoodHitsMap.at(pMCParticle);
	//	//	std::cout << "  " << std::endl;
	//	std::cout << "Primary MCParticle " << pMCParticle << std::endl;
	//	std::cout << "  - PDG : " << pMCParticle->GetParticleId() << std::endl;
	//std::cout << "  - NHits : " << caloHitList2.size() << std::endl; 
	Size = caloHitList2.size();
	mchitslist.push_back(std::make_pair(pMCParticle, Size));
	//	std::cout << "  " << std::endl;
      
      }



    
    for (const auto &ppfo : pfoVector)
      {
	mcpairslist.clear();
	
	const auto &caloHitList3 = PfosToGoodHitsMap.at(ppfo);
	//	std::cout << "Pfo " << ppfo << std::endl;
	//	std::cout << "  - PDG : " << ppfo->GetParticleId() << std::endl;
	//	std::cout << "  - NHits shared total: " << caloHitList3.size() << std::endl;  //this rewrites every time
	
	const auto iter = PfotoMCParticleMap.find(ppfo);
	if (iter == PfotoMCParticleMap.end()) {
	  std::cout << "Pfo not found." << std::endl;
	}
	else if ( (iter->second.size()) != 0) {
	  const auto hitsvector = iter->second;
	  for (int i = 0; i < hitsvector.size(); i++){
	    mcpairslist.push_back(std::make_pair(hitsvector[i].first, hitsvector[i].second.size()));
	  }
	}
     

	for (const Pfo *const pPPfo : *pPfoList) {
	  //  bool isClearCosmic2(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPPfo));
	  //METRIC
	  if (pPPfo == ppfo) {
	    //Get details about the pfo - need the list of hits in it
	    CaloHitList caloHitList;
	    LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_U, caloHitList);
	    LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_V, caloHitList);
	    LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_W, caloHitList);
	    unsigned int hitsinpfo = caloHitList.size();
	    float PurityTot = 0;
	    float SigTot = 0;
	    //bool isDownward = 1; //default as 1, will need to change to 0 (i.e. definetly a cosmic) //THINK HERE
	    //  if (hitsinpfo => 15) {
	      std::cout << "hits in pfo = " << hitsinpfo << std::endl;
	      for(it=mchitslist.begin(); it!=mchitslist.end(); ++it) {
		if (caloHitList3.size() != 0) {
		  std::pair<const pandora::MCParticle *, int> mchitsvalue = *it;
		  for(it2=mcpairslist.begin(); it2!=mcpairslist.end(); ++it2) {
		    std::pair<const pandora::MCParticle *, int> mcpairvalue = *it2;
		    if (mchitslist.size() == 1){
		      if (mchitsvalue.first == mcpairvalue.first && mcpairvalue.second != 0 && caloHitList3.size()==mcpairvalue.second) { //check mcvalue actually in this pfo	
			float puritytemp  = (float)caloHitList3.size()/(float)hitsinpfo;  
			float sigtemp = (float)mcpairvalue.second/(float)mchitsvalue.second;
			PurityTot = puritytemp;
			SigTot = SigTot + sigtemp;
			//	std::cout << "Purity temp = " << (float)caloHitList3.size() << "  over  " <<(float)hitsinpfo <<  std::endl;
			//	std::cout << "Sig temp = " << (float)mcpairvalue.second << "  over  " <<(float)mchitsvalue.second <<  std::endl;
			//----------------------------------------------------------------------------
			//	isDownward = mchitsvalue.first->GetMomentum().GetUnitVector().GetDotProduct(CartesianVector(0.f, 1.f, 0.f)) < 0;
    
		      }
		    }
		    else if (mchitsvalue.first == mcpairvalue.first && mcpairvalue.second != 0) {
		      const auto iterb = ParticletoPfoMap.find(mchitsvalue.first);   //look at particle to pfo map for a given mcparticle
		      if (iterb == ParticletoPfoMap.end()){
			std::cout << "MCParticle not found." << std::endl;
		      }
		      else if ( (iterb->second.size()) != 0) {      //if it is found with something in the vector
			const auto hitsvectorb = iterb->second;      //pull out the vector
			for (int i = 0; i < hitsvectorb.size(); i++){    //for each element of the vector
			  if (hitsvectorb[i].first == ppfo && hitsvectorb[i].second.size() != 0 ) {  //if pfo in vector matches current pfo + some number of shared hits is not 0
			    float puritytemp  = (float)caloHitList3.size()/(float)hitsinpfo;   
			    float sigtemp = (float)hitsvectorb[i].second.size()/(float)mchitsvalue.second;
			    PurityTot = puritytemp;
			    SigTot = SigTot + sigtemp;
			    // std::cout << "Purity temp = " << (float)caloHitList3.size() << "  over  " <<(float)hitsinpfo <<  std::endl;
			    // std::cout << "Sig temp = " << (float)hitsvectorb[i].second.size() << "  over  " <<(float)mchitsvalue.second <<  std::endl;
			    //  isDownward = mchitsvalue.first->GetMomentum().GetUnitVector().GetDotProduct(CartesianVector(0.f, 1.f, 0.f)) < 0;
			  }
			}

		      }
		    }
		  }
		}
	      }
	    
	      // // std::cout << "Purity = " << PurityTot << std::endl;
		// std::cout << "Sig = " << SigTot << std::endl;
	      // std::cout << "is down? " << isDownward << std::endl;

	      if (PurityTot < 0.05 && SigTot < 0.1){
		//	std::cout << "Clear Cosmic" << std::endl;

		clearcosmiccount = clearcosmiccount + 1;

		if (ppfo->GetParticleId() == 13) {
		  //	  std::cout << "Tagged as a CR by overall mechanism" << std::endl;
		  correctcr = correctcr + 1;
		}
		if (ppfo->GetParticleId()  == 12 || ppfo->GetParticleId()  == 14) {
		  //	  std::cout << "tagged as neutrino by overall mechanism" << std::endl;
		  incorrectcr = incorrectcr + 1;
		}
	      }
	      else if (PurityTot > 0.95 && SigTot > 0.1) {
		//	std::cout << "Clear Neutrino" << std::endl;

		clearneutrinocount = clearneutrinocount + 1;

		if (ppfo->GetParticleId() == 13) {
		  // std::cout << "Tagged as a CR by overall mechanism" << std::endl;
		  incorrectneutrino = incorrectneutrino + 1;
		}
		if (ppfo->GetParticleId()  == 12 || ppfo->GetParticleId()  == 14) {
		  // std::cout << "neutrino by overall mechanism" << std::endl;
		  correctneutrino = correctneutrino + 1;
		}
	      }
	      else{
		//	std::cout << "Unclear" << std::endl;

		if (ppfo->GetParticleId() == 13) {
		  //  std::cout << "Tagged as a CR by overall mechanism" << std::endl;
		  // incorrectneutrino = incorrectneutrino + 1;
		}
		if (ppfo->GetParticleId()  == 12 || ppfo->GetParticleId()  == 14) {
		  // std::cout << "neutrino by overall mechanism" << std::endl;
		  // correctneutrino = correctneutrino + 1;
		}
		    
	      }
	      // std::cout << "                                  " << std::endl;
	      
	      //  }
	      // else if (hitsinpfo==0){
	      if (hitsinpfo==0){
		std::cout << "This has been tagged as a neutrino." << std::endl;
		taggedneutrino = taggedneutrino + 1;
	      
	      }
	      // else{
	      // std::cout << "Number of hits too low to consider" << std::endl;
	      //  continue;
	      // }
		
	  }
   
	}
	
	std::cout << "                                  " << std::endl;
      }
  
    std::cout << "    " << std::endl;
    std::cout << "clear cosmic count " << clearcosmiccount << std::endl;
    std::cout << "clear neutrino count " << clearneutrinocount << std::endl;
    std::cout << "correct cr " << correctcr << std::endl;
    std::cout << "incorrect cr " << incorrectcr << std::endl;
    std::cout << "tagged nu " << taggedneutrino << std::endl;
    std::cout << "incorrect nu " << incorrectneutrino << std::endl;
    std::cout << "    " << std::endl;
    if ((incorrectneutrino == 0 && taggedneutrino == 0 && incorrectcr ==0)||(incorrectneutrino == 0 && taggedneutrino == 1 && incorrectcr == 0)) {
      std::cout << "Probably tagged correct neutrino/identified lack of neutrinos." << std::endl;
      wasright = 1;
    }
    else {
      std::cout << "Probably didn't find the right neutrino..." << std::endl;
      wasright = 0;
    } 
    std::cout << "    " << std::endl;
    std::cout << "was right? " << wasright << std::endl;
    std::cout << "    " << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreec", "wasright", wasright));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreec"));

    //Time to test the efficiency!!! 
    //need fraction of cosmic tagged/total cosmic as found by purity/sig
    std::cout << "   " << std::endl;
    float cosmiceff = (float)correctcr/(float)clearcosmiccount;
    std::cout << "Tagged clear cosmic = " << cosmiceff*100 << "%" << std::endl;
    float cosmicmis = (float)incorrectcr/(float)clearcosmiccount;
    std::cout << "clear cosmic identified as neutrino = " << cosmicmis*100 << "%" << std::endl;
    
    float neutrinoeff = (float)correctneutrino/(float)clearneutrinocount;
    std::cout << "protected clear neutrinos = " << neutrinoeff*100 << "%" << std::endl;
    float neutrinomis = (float)incorrectneutrino/(float)clearneutrinocount;
    std::cout << "clear neutrino identified as cosmic = " << neutrinomis*100 << "%" << std::endl;
    
  

    return STATUS_CODE_SUCCESS;
  }


  StatusCode PostMasterAnalysis::ReadSettings(const TiXmlHandle xmlHandle)
  {
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
													 "pfoListName", m_pfoListName));

    
    return STATUS_CODE_SUCCESS;
  }

  
}
