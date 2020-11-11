/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Api/PandoraContentApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

using namespace pandora;

namespace lar_content
{

MasterAlgorithm::MasterAlgorithm() :
    m_workerInstancesInitialized(false),
    m_larCaloHitVersion(1),
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(false),
    m_visualizeOverallRecoStatus(false),
    m_shouldRemoveOutOfTimeHits(true),
    m_pSlicingWorkerInstance(nullptr),
    m_pSliceNuWorkerInstance(nullptr),
    m_pSliceCRWorkerInstance(nullptr),
    m_fullWidthCRWorkerWireGaps(true),
    m_passMCParticlesToWorkerInstances(false),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_inTimeMaxX0(1.f)
{

}

  MasterAlgorithm::~MasterAlgorithm() {
    try
      {            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree", "output.root", "UPDATE"));
	//  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree2", "output2.root", "UPDATE"));
	//   PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree3", "output3.root", "UPDATE"));
	// PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree4", "output4.root", "UPDATE"));
	// PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttree5", "output5.root", "UPDATE"));
      }
    catch (const StatusCodeException &)
      {
	std::cout << " Unable to write tree  to file " << std::endl;
      }
  }

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::ShiftPfoHierarchy(const ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const
{
    if (!pParentPfo->GetParentPfoList().empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PfoToLArTPCMap::const_iterator larTPCIter(pfoToLArTPCMap.find(pParentPfo));

    if (pfoToLArTPCMap.end() == larTPCIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);


    PfoList pfoList;
    LArPfoHelper::GetAllDownstreamPfos(pParentPfo, pfoList);

    if (m_visualizeOverallRecoStatus)
    {
        std::cout << "ShiftPfoHierarchy: x0 " << x0 << std::endl;
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "BeforeShiftCRPfos", GREEN));
    }

    for (const ParticleFlowObject *const pDaughterPfo : pfoList)
    {
        CaloHitList caloHitList;
        for (const Cluster *const pCluster : pDaughterPfo->GetClusterList())
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        }

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            PandoraContentApi::CaloHit::Metadata metadata;
            metadata.m_x0 = x0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
        }

        for (const Vertex *const pVertex : pDaughterPfo->GetVertexList())
        {
            PandoraContentApi::Vertex::Metadata metadata;
            metadata.m_x0 = x0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::AlterMetadata(*this, pVertex, metadata));
        }
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "AfterShiftCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::StitchPfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete,
    PfoToLArTPCMap &pfoToLArTPCMap) const
{
    if (pPfoToEnlarge == pPfoToDelete)
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN Remove pfos from pfo to lar tpc map here, avoiding problems if stitching across multiple tpcs.
    pfoToLArTPCMap.erase(pPfoToEnlarge);
    pfoToLArTPCMap.erase(pPfoToDelete);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, m_recreatedPfoListName));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));

    for (const  Vertex *const pDaughterVertex : daughterVertices)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, m_recreatedVertexListName));

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(PfoMopUpBaseAlgorithm::GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
                m_recreatedClusterListName, m_recreatedClusterListName));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Run()
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Reset());

    if (!m_workerInstancesInitialized)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->InitializeWorkerInstances());
    

    if (m_passMCParticlesToWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CopyMCParticles());


    PfoToFloatMap stitchedPfosToX0Map;
    VolumeIdToHitListMap volumeIdToHitListMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVolumeIdToHitListMap(volumeIdToHitListMap));


    if (m_shouldRunAllHitsCosmicReco)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayReconstruction(volumeIdToHitListMap));


        PfoToLArTPCMap pfoToLArTPCMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RecreateCosmicRayPfos(pfoToLArTPCMap));


        if (m_shouldRunStitching)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->StitchCosmicRayPfos(pfoToLArTPCMap, stitchedPfosToX0Map));

    }

    if (m_shouldRunCosmicHitRemoval)
    {
        PfoList clearCosmicRayPfos, ambiguousPfos;
	bool directioncosmic;
	float downprobability = -1;
	//	float deltachi2;
	//	float deltachi2alone;
	//	float minchi2perhit;
	//	bool incomplete;
	// PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos, directioncosmic, downprobability, deltachi2, deltachi2alone, minchi2perhit, incomplete));
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(stitchedPfosToX0Map, clearCosmicRayPfos, ambiguousPfos, directioncosmic, downprobability));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayHitRemoval(ambiguousPfos));

    }

    SliceVector sliceVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSlicing(volumeIdToHitListMap, sliceVector));

    if (m_shouldRunNeutrinoRecoOption || m_shouldRunCosmicRecoOption)
    {
        SliceHypotheses nuSliceHypotheses, crSliceHypotheses;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSliceReconstruction(sliceVector, nuSliceHypotheses, crSliceHypotheses));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SelectBestSliceHypotheses(nuSliceHypotheses, crSliceHypotheses, sliceVector));


    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::InitializeWorkerInstances()
{
    // ATTN Used to be in the regular Initialize callback, but detector gap list cannot be extracted in client app before the first event
    if (m_workerInstancesInitialized)
        return STATUS_CODE_ALREADY_INITIALIZED;

    try
    {
        const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
        const DetectorGapList &gapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

        for (const LArTPCMap::value_type &mapEntry : larTPCMap)
        {
            const unsigned int volumeId(mapEntry.second->GetLArTPCVolumeId());
            m_crWorkerInstances.push_back(this->CreateWorkerInstance(*(mapEntry.second), gapList, m_crSettingsFile, "CRWorkerInstance" + std::to_string(volumeId)));
        }

        if (m_shouldRunSlicing)
            m_pSlicingWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_slicingSettingsFile, "SlicingWorker");

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile, "SliceNuWorker");

        if (m_shouldRunCosmicRecoOption)
            m_pSliceCRWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile, "SliceCRWorker");
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "MasterAlgorithm: Exception during initialization of worker instances " << statusCodeException.ToString() << std::endl;
        return statusCodeException.GetStatusCode();
    }

    m_workerInstancesInitialized = true;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::CopyMCParticles() const
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    PandoraInstanceList pandoraWorkerInstances(m_crWorkerInstances);
    if (m_pSlicingWorkerInstance) pandoraWorkerInstances.push_back(m_pSlicingWorkerInstance);
    if (m_pSliceNuWorkerInstance) pandoraWorkerInstances.push_back(m_pSliceNuWorkerInstance);
    if (m_pSliceCRWorkerInstance) pandoraWorkerInstances.push_back(m_pSliceCRWorkerInstance);

    LArMCParticleFactory mcParticleFactory;

    for (const Pandora *const pPandoraWorker : pandoraWorkerInstances)
    {
        for (const MCParticle *const pMCParticle : *pMCParticleList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pPandoraWorker, pMCParticle, &mcParticleFactory));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const unsigned int nLArTPCs(larTPCMap.size());

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit*>(pCaloHit));

        if (!pLArCaloHit && (1 != nLArTPCs))
            return STATUS_CODE_INVALID_PARAMETER;

        const unsigned int volumeId(pLArCaloHit ? pLArCaloHit->GetLArTPCVolumeId() : 0);
        const LArTPC *const pLArTPC(larTPCMap.at(volumeId));

        LArTPCHitList &larTPCHitList(volumeIdToHitListMap[volumeId]);
        larTPCHitList.m_allHitList.push_back(pCaloHit);

        if (((pCaloHit->GetPositionVector().GetX() >= (pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX())) &&
            (pCaloHit->GetPositionVector().GetX() <= (pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX()))))
        {
            larTPCHitList.m_truncatedHitList.push_back(pCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    unsigned int workerCounter(0);

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());
        VolumeIdToHitListMap::const_iterator iter(volumeIdToHitListMap.find(larTPC.GetLArTPCVolumeId()));

        if (volumeIdToHitListMap.end() == iter)
            continue;

        for (const CaloHit *const pCaloHit : iter->second.m_allHitList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pCRWorker, pCaloHit));

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance " << ++workerCounter << " of " << m_crWorkerInstances.size() << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pCRWorker));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RecreateCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const PfoList *pCRPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pCRWorker, pCRPfos));

        PfoList newPfoList;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(*pCRPfos, newPfoList));

        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());

        for (const Pfo *const pNewPfo : newPfoList)
            pfoToLArTPCMap[pNewPfo] = &larTPC;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::StitchCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "RecreatedCRPfos", GREEN));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    for (StitchingBaseTool *const pStitchingTool : m_stitchingToolVector)
        pStitchingTool->Run(this, pRecreatedCRPfos, pfoToLArTPCMap, stitchedPfosToX0Map);

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "AfterStitchingCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// StatusCode MasterAlgorithm::TagCosmicRayPfos(const PfoToFloatMap &stitchedPfosToX0Map, PfoList &clearCosmicRayPfos, PfoList &ambiguousPfos, bool &directioncosmic, float &downprobability, float &deltachi2, float &deltachi2alone, float &minchi2perhit, bool &incomplete) const
 StatusCode MasterAlgorithm::TagCosmicRayPfos(const PfoToFloatMap &stitchedPfosToX0Map, PfoList &clearCosmicRayPfos, PfoList &ambiguousPfos, bool &directioncosmic, float &downprobability) const
  {
    //UNDO FOR TDT
    //-------------------------------------------------------down
    //METRIC
    //int Size = 0; //mc stuff
    std::list<std::pair<const pandora::MCParticle*, int>> mchitslist;
    std::list<std::pair<const pandora::MCParticle*, int>> mcpairslist;
    MCParticleList mcparticlelist;
    std::list<std::pair<const pandora::MCParticle*, int>>::iterator it;
    std::list<std::pair<const pandora::MCParticle*, int>>::iterator it2;
    CaloHitList totalcalohits;
    //PfoList totalpfolist;
    // int clearcosmiccount = 0; //mc stuff
    // int clearneutrinocount = 0;
    // int correctcr = 0;
    // int incorrectcr = 0;
    //int correctneutrino = 0;
    //int incorrectneutrino = 0;
    //bool isClearCosmic = false;
    std::vector<float> viewsizelist;
 
    //------------------------------------------------------------------up
    //---------------------------------------------------------------------------------------------------------
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    PfoList nonStitchedParentCosmicRayPfos;
    for (const Pfo *const pPfo : *pRecreatedCRPfos)
      {
	if (!pPfo->GetParentPfoList().empty())
	  continue;

	PfoToFloatMap::const_iterator pfoToX0Iter = stitchedPfosToX0Map.find(pPfo);
	const float x0Shift((pfoToX0Iter != stitchedPfosToX0Map.end()) ? pfoToX0Iter->second : 0.f);
	PfoList &targetList((std::fabs(x0Shift) > m_inTimeMaxX0) ? clearCosmicRayPfos : nonStitchedParentCosmicRayPfos);
	targetList.push_back(pPfo);
      }

    for (CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool : m_cosmicRayTaggingToolVector)
      pCosmicRayTaggingTool->FindAmbiguousPfos(nonStitchedParentCosmicRayPfos, ambiguousPfos, this);

    

    //UNDO
    //-------------------------------------------------------------down
    // std::unordered_map<const Pfo*, float> pfotoprobabilitymap;
    /*
    std::unordered_map<const Pfo*, float> pfotodeltachi2map;
    std::unordered_map<const Pfo*, float> pfotodeltachi2alonemap;
    std::unordered_map<const Pfo*, float> pfotominchi2perhitmap;
    */
    // std::unordered_map<const Pfo*, float> truecosmicpfotoprobabilitymap;
    // std::unordered_map<const Pfo*, float> truecosmicpfotodeltachi2map;
    // std::unordered_map<const Pfo*, float> truecosmicpfotodeltachi2alonemap;
    //std::unordered_map<const Pfo*, bool>  truecosmicpfotoisdownmap;
    // std::unordered_map<const Pfo*, float>  truecosmicpfotominchi2perhitmap;
    //  std::unordered_map<const Pfo*, float> trueneutrinopfotoprobabilitymap;
    // std::unordered_map<const Pfo*, float> trueneutrinopfotodeltachi2map;
    // std::unordered_map<const Pfo*, float> trueneutrinopfotodeltachi2alonemap;
    // std::unordered_map<const Pfo*, bool>  trueneutrinopfotoisdownmap;
    // std::unordered_map<const Pfo*, float>  trueneutrinopfotominchi2perhitmap;
    //---------------------------------------------------------------------up

    /*
    for (const Pfo *const pPfo : nonStitchedParentCosmicRayPfos)
    {
        const bool isClearCosmic(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo));

        if (isClearCosmic)
            clearCosmicRayPfos.push_back(pPfo);
    }
    */

    //UNDO
    //---------------------------------------------------------------------------------down


    for (const Pfo *const pPPPfo : nonStitchedParentCosmicRayPfos) {

      bool isClearCosmic(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPPPfo));
      directioncosmic = false;
      downprobability = -1.0;
      if(!isClearCosmic){  //if it's not already clearly a cosmic, see if it's above 80% probability of being one
	//std::cout << "     " << std::endl;
	//std::cout << "CLEAR FINAL STATE? " << LArPfoHelper::IsFinalState(pPPPfo) << std::endl;
	/*
	  for (TrackDirectionBaseTool *const pTrackDirectionTool : m_trackDirectionToolVector){
	  try{
	  // pTrackDirectionTool->FindDirections(pPPPfo, directioncosmic, downprobability, deltachi2, deltachi2alone, minchi2perhit, incomplete, this);  //probability
	  pTrackDirectionTool->FindDirections(pPPPfo, directioncosmic, downprobability, this);  //probability 
	  }
	  catch(...){
	  std::cout << "Runtime error in TrackDirectionTool!" << std::endl;
	  }
	*/
	// if (incomplete == false) {
	// pfotoprobabilitymap.insert({pPPPfo,downprobability});
	/*
	  pfotodeltachi2map.insert({pPPPfo,deltachi2});
	  pfotodeltachi2alonemap.insert({pPPPfo,deltachi2alone});
	  pfotominchi2perhitmap.insert({pPPPfo, minchi2perhit});
	*/
	// std::cout << "bottom of the try block" << std::endl;
	// }

//----------CUTS-----------------------------------------------------------------------------------
	/*
	CaloHitList totalcalohitsN;
	CaloHitList totalcalohitsW;
	LArPfoHelper::GetCaloHits(pPPPfo, TPC_VIEW_U, totalcalohitsN);
	LArPfoHelper::GetCaloHits(pPPPfo, TPC_VIEW_V, totalcalohitsN);
	LArPfoHelper::GetCaloHits(pPPPfo, TPC_VIEW_W, totalcalohitsN);
	LArPfoHelper::GetCaloHits(pPPPfo, TPC_VIEW_W, totalcalohitsW);
	std::cout << "total = " << totalcalohitsN.size() << std::endl;
	std::cout << "total W = " << totalcalohitsW.size() << std::endl;


	ClusterList clusterList = pPPPfo->GetClusterList();
	std::cout << "cluster size = " << clusterList.size() << std::endl;
	const Cluster* pCluster = clusterList.back();

	int w = 0;
	for (auto c :clusterList) {
	  if (w == 2) {
	    pCluster = c;
	  }
	  w++;
	}
	    
	float length = LArClusterHelper::GetLength(pCluster);
	std::cout << "length = " << length << std::endl;

	
	if(totalcalohitsN.size() < 15) {
	  directioncosmic = true;
	  downprobability = -2.0;
	}	    
	else if (totalcalohitsW.size() < 4 && totalcalohitsW.size() > 0 ) {
	  directioncosmic = true;
	  downprobability = -3.0;
	}	   
	else if(length < 1.0 && totalcalohitsW.size() != 0) {
	  directioncosmic = true;
	  downprobability = -4.0;
	      
	}
	
	    
	  
	std::cout << "final directioncosmic : " << directioncosmic << std::endl;
	std::cout << "down probability : " << downprobability << std::endl;

	if(directioncosmic == true) {
	  isClearCosmic = true;               //if it is, change the bool
	}

      }
	*/
      }
//-----------------CUTS----------------------------------------------------------------------
      if (isClearCosmic) {
	clearCosmicRayPfos.push_back(pPPPfo);
	//	ambiguousPfos.remove(pPPPfo);    // could do this later, if it's already been found as a pfo it won't be again
      }
    }

    //}
  //----------------------------------------------------------------------------------------------------up

  //  if (incomplete == false) {
      //UNDO
      //-------------------------------------------------------------------------------------------------down
    /*  Top of MC things
    std::cout << "  " << std::endl;
    std::cout << "*pRecreatedCRPfos.size " << pRecreatedCRPfos->size() << std::endl;
      for (const Pfo *const pPPPPfo : *pRecreatedCRPfos ) {
	try{
	  const MCParticle *pMCParticle = LArMCParticleHelper::GetMainMCParticle(pPPPPfo);
	  mcparticlelist.push_back(pMCParticle);
	  //totalpfolist.push_back(pPPPfo);
	  LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_U, totalcalohits);
	  LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_V, totalcalohits);
	  LArPfoHelper::GetCaloHits(pPPPPfo, TPC_VIEW_W, totalcalohits);
	  std::cout << "MC particle found" << std::endl;
	}
	catch(...){
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

      std::cout << "mcMCParticlesToGoodHitsMap,size() " << mcMCParticlesToGoodHitsMap.size() << std::endl;


  
      //get the pfos to hits map
      LArMCParticleHelper::PfoContributionMap PfosToGoodHitsMap;
      LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pRecreatedCRPfos, {mcMCParticlesToGoodHitsMap}, PfosToGoodHitsMap);  //changed

 

      //find hits shared between Pfos and MCParticles
      LArMCParticleHelper::PfoToMCParticleHitSharingMap PfotoMCParticleMap;
      LArMCParticleHelper::MCParticleToPfoHitSharingMap ParticletoPfoMap;
      LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(PfosToGoodHitsMap, {mcMCParticlesToGoodHitsMap}, PfotoMCParticleMap, ParticletoPfoMap);


    
      MCParticleVector mcParticleVector;
      LArMonitoringHelper::GetOrderedMCParticleVector({mcMCParticlesToGoodHitsMap}, mcParticleVector);
      PfoVector pfoVector;
      LArMonitoringHelper::GetOrderedPfoVector({PfosToGoodHitsMap}, pfoVector);
      


      unsigned int nMatches = std::numeric_limits<unsigned int>::max();
      std::cout << "nMatches = " << nMatches << std::endl;
      LArMonitoringHelper::PrintMatchingTable(pfoVector, mcParticleVector, ParticletoPfoMap, nMatches);
      std::cout << "nMatches = " << nMatches << std::endl;


    
      for (const auto &pMCParticle : mcParticleVector)
	{
	  const auto &caloHitList2 = mcMCParticlesToGoodHitsMap.at(pMCParticle);
	  // std::cout << "  " << std::endl;
	  // std::cout << "Primary MCParticle " << pMCParticle << std::endl;
	  // std::cout << "  - PDG : " << pMCParticle->GetParticleId() << std::endl;
	  // std::cout << "  - NHits : " << caloHitList2.size() << std::endl; 
	  Size = caloHitList2.size();
	  mchitslist.push_back(std::make_pair(pMCParticle, Size));
	  // std::cout << "  " << std::endl;
      
	}

    
      for (const auto &ppfo : pfoVector)
	{
	  mcpairslist.clear();
	
	  const auto &caloHitList3 = PfosToGoodHitsMap.at(ppfo);
	  // std::cout << "Pfo " << ppfo << std::endl;
	  // std::cout << "  - PDG : " << ppfo->GetParticleId() << std::endl;
	  // std::cout << "  - NHits shared total: " << caloHitList3.size() << std::endl;
	  //  std::cout << "  - is final state?: " << LArPfoHelper::IsFinalState(ppfo)  << std::endl;
	
	  const auto iter = PfotoMCParticleMap.find(ppfo);
	  // if (iter == PfotoMCParticleMap.end()) {
	  //   std::cout << "Pfo not found." << std::endl;
	  // }
	  if ( (iter->second.size()) != 0) {
	    const auto hitsvector = iter->second;
	    for (int i = 0; i < hitsvector.size(); i++){
	      mcpairslist.push_back(std::make_pair(hitsvector[i].first, hitsvector[i].second.size()));
	    }
	  }

     
      
	  for (const Pfo *const pPPfo : *pRecreatedCRPfos) {
	    if(LArPfoHelper::IsFinalState(pPPfo)) {
	      //bool isNotClearCosmic2(clearCosmicRayPfos.end() == std::find(clearCosmicRayPfos.begin(), clearCosmicRayPfos.end(), pPPfo));
	      bool isClearCosmic2(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPPfo));
	      //METRIC

	      // for (auto& y: pfoToIsLikelyCRMuonMap ) {
	      //  if (x.first == ppfo){
	      if (pPPfo == ppfo) {
		//Get details about the pfo - need the list of hits in it
		CaloHitList caloHitList;
		CaloHitList caloHitListU;
		CaloHitList caloHitListV;
		CaloHitList caloHitListW;
		LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_U, caloHitList);
		LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_V, caloHitList);
		LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_W, caloHitList);
		LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_U, caloHitListU);
		LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_V, caloHitListV);
		LArPfoHelper::GetCaloHits(pPPfo, TPC_VIEW_W, caloHitListW);
      
	        float hitsinpfo = caloHitList.size();
    */ //Bottom of MC 1
		/*
		ClusterList clusterList = ppfo->GetClusterList();
		
		const Cluster* pCluster = clusterList.back();
		int w = 0;
		for (auto c :clusterList) {
		  if (w == 2) {
		    pCluster = c;
		  }
		  w++;
		}
		*/
		//	float length = LArClusterHelper::GetLength(pCluster);
    /* MC stuff 2
		if (caloHitListU.size() >= 5) {
		  viewsizelist.push_back(caloHitListU.size());
		}
		if (caloHitListV.size() >= 5) {
		  viewsizelist.push_back(caloHitListV.size());
		}
		if (caloHitListW.size() >= 5) {
		  viewsizelist.push_back(caloHitListW.size());
		}

		

		float PurityTot = 0;
		float SigTot = 0;
		//	bool isDownward = 1; //default as 1, will need to change to 0 (i.e. definetly a neutrino)
		if (hitsinpfo >= 15) {
		  if (viewsizelist.size() >= 2) {
		    //  std::cout << "hits in pfo = " << hitsinpfo << std::endl;
		    //  std::cout << "hits in pfo W = " << caloHitListW.size() << std::endl;
		    //	  std::cout << "length in W = " << length << std::endl;
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
			    // std::cout << "Purity temp = " << (float)caloHitList3.size() << "  over  " <<(float)hitsinpfo <<  std::endl;
			    //  std::cout << "Sig temp = " << (float)mcpairvalue.second << "  over  " <<(float)mchitsvalue.second <<  std::endl;
			    //----------------------------------------------------------------------------
			    // isDownward = mchitsvalue.first->GetMomentum().GetUnitVector().GetDotProduct(CartesianVector(0.f, 1.f, 0.f)) < 0;
			    //std::cout << "is down? " << isDownward << std::endl;
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
				//	std::cout << "Purity temp = " << (float)caloHitList3.size() << "  over  " <<(float)hitsinpfo <<  std::endl;
				//	std::cout << "Sig temp = " << (float)hitsvectorb[i].second.size() << "  over  " <<(float)mchitsvalue.second <<  std::endl;
				//  isDownward = mchitsvalue.first->GetMomentum().GetUnitVector().GetDotProduct(CartesianVector(0.f, 1.f, 0.f)) < 0;
				//std::cout << "is down? " << isDownward << std::endl;
			      }
			    }

			  }
			}
		      }
		    }
		  }
	    
		  // std::cout << "Purity = " << PurityTot << std::endl;
		  //  std::cout << "Sig = " << SigTot << std::endl;
		  //  std::cout << "is down? " << isDownward << std::endl;

		  if (PurityTot < 0.05 && SigTot < 0.1){
		    //std::cout << "Clear Cosmic" << std::endl;
		    clearcosmiccount = clearcosmiccount + 1;
	        
		    if (isClearCosmic2 == 1) {
		      //  std::cout << "Tagged as a CR" << std::endl;
		      correctcr = correctcr + 1;
		    }
		    if (isClearCosmic2 == 0) {
		      // std::cout << "potential neutrino" << std::endl;
		      incorrectcr = incorrectcr + 1;
		    }
		  }
		  else if (PurityTot > 0.95 && SigTot > 0.1) {
		    // std::cout << "Clear Neutrino" << std::endl;

		    clearneutrinocount = clearneutrinocount + 1;
		    if (isClearCosmic2 == 1) {
		      //  std::cout << "Tagged as a CR" << std::endl;
		      incorrectneutrino = incorrectneutrino + 1;
		    }
		    if (isClearCosmic2 == 0) {
		      //  std::cout << "potential neutrino" << std::endl;
		      correctneutrino = correctneutrino + 1;
		    }
		  }
		  else{
		    //  std::cout << "Unclear" << std::endl;
		  }
		  //std::cout << "                                  " << std::endl;
	      
		  }
		}
		else{
		  // std::cout << "Number of hits too low to consider" << std::endl;
		  continue;
		}
	      }
   
	    }
	  
	  }
	}
  

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
    

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "cosmiceff", cosmiceff));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "cosmicmis", cosmicmis));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "neutrinoeff", neutrinoeff));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttree", "neutrinomis", neutrinomis));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttree"));
    
    */ //Mc stuff 2
    
    // PANDORA_MONITORING_API(ScanTree(this->GetPandora(), "ttree"))

    //for (const Pfo *const pPfo : clearCosmicRayPfos) {
    // bool found(ambiguousPfos.end() != std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo));
    // if (found == true) {
    //	ambiguousPfos.remove(pPfo);
	//}
    // }
    //  std::cout << "SIZE END = " << ambiguousPfos.size() << std::endl;
    //---------------------------------------------------------------------------------------------------------------up
    //--------------------------------------------------------------------------------------------------
    for (const Pfo *const pPfo : *pRecreatedCRPfos) ///repeat the process for a different set
      {
        bool isClearCosmic(ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo));
	//	if(!isClearCosmic){       //BACK IN LATER - don't think need now has have removed all found clear cosmics from ambiguousPfos
	//	  for (TrackDirectionBaseTool *const pTrackDirectionTool : m_trackDirectionToolVector){
	//	    pTrackDirectionTool->FindDirections(pPfo, directioncosmic, this);
	//std::cout << "final directioncosmic : " << directioncosmic << std::endl;

	//    if(directioncosmic == true) {
	//      isClearCosmic = true;
	//	    }

	//	  }
	//	}
	
        PandoraContentApi::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["IsClearCosmic"] = (isClearCosmic ? 1.f : 0.f);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
      }

    if (m_visualizeOverallRecoStatus)
      {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &clearCosmicRayPfos, "ClearCRPfos", RED));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &ambiguousPfos, "AmbiguousCRPfos", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
      }

    return STATUS_CODE_SUCCESS;
  }

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayHitRemoval(const PfoList &ambiguousPfos) const
{
    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(ambiguousPfos, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        const ClusterList clusterList(pPfoToDelete->GetClusterList());
        const VertexList vertexList(pPfoToDelete->GetVertexList());

        // ATTN: If an ambiguous pfo has been stitched, reset the calo hit positions in preparation for subsequent algorithm chains
        if (LArStitchingHelper::HasPfoBeenStitched(pPfoToDelete))
        {
            CaloHitList caloHitList2D;
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfoToDelete, TPC_VIEW_W, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfoToDelete, TPC_VIEW_W, caloHitList2D);

            for (const CaloHit *const pCaloHit : caloHitList2D)
            {
                PandoraContentApi::CaloHit::Metadata metadata;
                metadata.m_x0 = 0.f;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
            }
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &vertexList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const
{
    for (const VolumeIdToHitListMap::value_type &mapEntry : volumeIdToHitListMap)
    {
        for (const CaloHit *const pCaloHit : (m_shouldRemoveOutOfTimeHits ? mapEntry.second.m_truncatedHitList : mapEntry.second.m_allHitList))
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            const HitType hitType(pCaloHit->GetHitType());
            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                continue;

            if (m_shouldRunSlicing)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSlicingWorkerInstance, pCaloHit));
            }
            else
            {
                if (sliceVector.empty()) sliceVector.push_back(CaloHitList());
                sliceVector.back().push_back(pCaloHit);
            }
        }
    }

    if (m_shouldRunSlicing)
    {
        if (m_printOverallRecoStatus)
            std::cout << "Running slicing worker instance" << std::endl;

        const PfoList *pSlicePfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSlicingWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSlicingWorkerInstance, pSlicePfos));

        if (m_visualizeOverallRecoStatus)
        {
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pSlicePfos, "OnePfoPerSlice", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        for (const Pfo *const pSlicePfo : *pSlicePfos)
        {
            sliceVector.push_back(CaloHitList());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_U, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_V, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_W, sliceVector.back());
        }
    }

    if (m_printOverallRecoStatus)
        std::cout << "Identified " << sliceVector.size() << " slice(s)" << std::endl;


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSliceReconstruction(SliceVector &sliceVector, SliceHypotheses &nuSliceHypotheses, SliceHypotheses &crSliceHypotheses) const
{
    SliceVector selectedSliceVector;
    if (m_shouldRunSlicing && !m_sliceSelectionToolVector.empty())
    {
        SliceVector inputSliceVector(sliceVector);
        for (SliceSelectionBaseTool *const pSliceSelectionTool : m_sliceSelectionToolVector)
        {
            pSliceSelectionTool->SelectSlices(this, inputSliceVector, selectedSliceVector);
            inputSliceVector = selectedSliceVector;
        }
    }
    else
    {
        selectedSliceVector = std::move(sliceVector);
    }

    unsigned int sliceCounter(0);

    for (const CaloHitList &sliceHits : selectedSliceVector)
    {
        for (const CaloHit *const pSliceCaloHit : sliceHits)
        {
            // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
            const CaloHit *const pCaloHitInMaster(m_shouldRunSlicing ? static_cast<const CaloHit*>(pSliceCaloHit->GetParentAddress()) : pSliceCaloHit);

            if (m_shouldRunNeutrinoRecoOption)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceNuWorkerInstance, pCaloHitInMaster));

            if (m_shouldRunCosmicRecoOption)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceCRWorkerInstance, pCaloHitInMaster));
        }

        if (m_shouldRunNeutrinoRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running nu worker instance for slice " << (sliceCounter + 1) << " of " << selectedSliceVector.size() << std::endl;

            const PfoList *pSliceNuPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceNuWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceNuWorkerInstance, pSliceNuPfos));
            nuSliceHypotheses.push_back(*pSliceNuPfos);

            for (const ParticleFlowObject *const pPfo : *pSliceNuPfos)
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["SliceIndex"] = sliceCounter;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
        }

        if (m_shouldRunCosmicRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running cr worker instance for slice " << (sliceCounter + 1) << " of " << selectedSliceVector.size() << std::endl;


            const PfoList *pSliceCRPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceCRWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceCRWorkerInstance, pSliceCRPfos));
            crSliceHypotheses.push_back(*pSliceCRPfos);

            for (const ParticleFlowObject *const pPfo : *pSliceCRPfos)
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["SliceIndex"] = sliceCounter;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
        }

        ++sliceCounter;
    }

    // ATTN: If we swapped these objects at the start, be sure to swap them back in case we ever want to use sliceVector
    // after this function
    if (!(m_shouldRunSlicing && !m_sliceSelectionToolVector.empty()))
        sliceVector = std::move(selectedSliceVector);

    if (m_shouldRunNeutrinoRecoOption && m_shouldRunCosmicRecoOption && (nuSliceHypotheses.size() != crSliceHypotheses.size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

  StatusCode MasterAlgorithm::SelectBestSliceHypotheses(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, const SliceVector &sliceVector) const
{
    if (m_printOverallRecoStatus)
        std::cout << "Select best slice hypotheses" << std::endl;
    
    
    PfoList selectedSlicePfos; ////giant pfo list
    
    
    PfoList selectedSlicePfosB;
    //  float downprobability = -1;
    bool directioncosmic;

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
       const PfoList &neutrinoPfoList(nuSliceHypotheses.at(sliceIndex));

       for (const Pfo *const pNeutrinoPfo : neutrinoPfoList)
	 {
	   PfoList daughterPfos = pNeutrinoPfo->GetDaughterPfoList();
	   for (const ParticleFlowObject *const pDaughterPfo : daughterPfos) {
	     //std::cout << " daughters " << pDaughterPfo  << std::endl;
	     selectedSlicePfosB.push_back(pDaughterPfo);
	   }

	 }

    }



    for (unsigned int sliceIndex = 0, nSlices = crSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
       const PfoList &cosmicPfoList(crSliceHypotheses.at(sliceIndex));

       for (const Pfo *const pcosmicPfo : cosmicPfoList)
        {
	  //std::cout << "crnut " << pcosmicPfo << std::endl;
	  selectedSlicePfosB.push_back(pcosmicPfo);

	}

    }

    //std::cout << " size   :   " <<  selectedSlicePfosB.size() << std::endl; 

    std::unordered_map<const Pfo*, float> pfotoprobabilitymap;
    PfoToFloatMap pfotoprobabilitymapb;
    //std::cout << "------------------Running over Slices---------------------" << std::endl;
    for (const Pfo *const pPPPfo :  selectedSlicePfosB) {
      float downprobability = -1;
      for (TrackDirectionBaseTool *const pTrackDirectionTool : m_trackDirectionToolVector){
	try{
	  pTrackDirectionTool->FindDirections(pPPPfo, directioncosmic, downprobability, this);  //probability
	  pfotoprobabilitymapb.insert({pPPPfo,downprobability});
	  // std::cout << " ^^ Pfo = " << pPPPfo << std::endl;
	}
	catch(...){
	  std::cout << "Runtime error in TrackDirectionTool!" << std::endl;
	}
      }
    }
    //  std::cout << "------------------Done Running over Slices---------------------" << std::endl;
    
    
    //probability from list
    //map of pfo to probability
    //map into select output pfos
    // PfoToFloatMap pfotoprobabilitymapb;

    if (m_shouldPerformSliceId)
    {
      // std::cout << "option 1 run SelectoutputPfos" << std::endl;
        for (SliceIdBaseTool *const pSliceIdTool : m_sliceIdToolVector)
	  pSliceIdTool->SelectOutputPfos(this, nuSliceHypotheses, crSliceHypotheses, selectedSlicePfos, pfotoprobabilitymapb, sliceVector);

    }
    else if (m_shouldRunNeutrinoRecoOption != m_shouldRunCosmicRecoOption)
    {
      // std::cout << "option 2 run call em crs" << std::endl;
        const SliceHypotheses &sliceHypotheses(m_shouldRunNeutrinoRecoOption ? nuSliceHypotheses : crSliceHypotheses);

        for (const PfoList &slice : sliceHypotheses)
            selectedSlicePfos.insert(selectedSlicePfos.end(), slice.begin(), slice.end());
    }

 
    PfoList newSlicePfoList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(selectedSlicePfos, newSlicePfoList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Reset()
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pCRWorker));

    if (m_pSlicingWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSlicingWorkerInstance));

    if (m_pSliceNuWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceNuWorkerInstance));

    if (m_pSliceCRWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceCRWorkerInstance));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const CaloHit *const pCaloHit) const
{
    const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit*>(pCaloHit)};
    if (pLArCaloHit == nullptr)
    {
        std::cout << "MasterAlgorithm: Could not cast CaloHit to LArCaloHit" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    LArCaloHitParameters parameters;
    parameters.m_positionVector = pCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pCaloHit->GetCellSize1();
    parameters.m_cellThickness = pCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit->GetTime();
    parameters.m_inputEnergy = pCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit->IsDigital();
    parameters.m_hitType = pCaloHit->GetHitType();
    parameters.m_hitRegion = pCaloHit->GetHitRegion();
    parameters.m_layer = pCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit->IsInOuterSamplingLayer();
    // ATTN Parent of calo hit in worker is corresponding calo hit in master
    parameters.m_pParentAddress = static_cast<const void*>(pCaloHit);
    parameters.m_larTPCVolumeId = pLArCaloHit->GetLArTPCVolumeId();
    parameters.m_daughterVolumeId = (m_larCaloHitVersion > 1) ? pLArCaloHit->GetDaughterVolumeId() : 0;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, parameters, m_larCaloHitFactory));

    if (m_passMCParticlesToWorkerInstances)
    {
        MCParticleVector mcParticleVector;
        for (const auto &weightMapEntry : pLArCaloHit->GetMCParticleWeightMap())
            mcParticleVector.push_back(weightMapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(
                *pPandora, pLArCaloHit, pMCParticle, pLArCaloHit->GetMCParticleWeightMap().at(pMCParticle)));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const MCParticle *const pMCParticle, const LArMCParticleFactory *const pMCParticleFactory) const
{
    LArMCParticleParameters parameters;
    const LArMCParticle *const pLArMCParticle = dynamic_cast<const LArMCParticle*>(pMCParticle);

    if (!pLArMCParticle)
    {
        std::cout << "MasterAlgorithm::Copy - Expect to pass only LArMCParticles to Pandora worker instances." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    parameters.m_nuanceCode = pLArMCParticle->GetNuanceCode();
    parameters.m_energy = pMCParticle->GetEnergy();
    parameters.m_momentum = pMCParticle->GetMomentum();
    parameters.m_vertex = pMCParticle->GetVertex();
    parameters.m_endpoint = pMCParticle->GetEndpoint();
    parameters.m_particleId = pMCParticle->GetParticleId();
    parameters.m_mcParticleType = pMCParticle->GetMCParticleType();
    // ATTN Parent of mc particle in worker is corresponding mc particle in master
    parameters.m_pParentAddress = static_cast<const void*>(pMCParticle);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, parameters, *pMCParticleFactory));

    for (const MCParticle *const pDaughterMCParticle : pMCParticle->GetDaughterList())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pMCParticle, pDaughterMCParticle));

    for (const MCParticle *const pParentMCParticle : pMCParticle->GetParentList())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora, pParentMCParticle, pMCParticle));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const PfoList &inputPfoList, PfoList &newPfoList) const
{
    if (inputPfoList.empty())
        return STATUS_CODE_SUCCESS;

    // TODO if no pfo in input list is primary - raise exception

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Pfo *const pPfo : inputPfoList)
    {
        if (pPfo->GetParentPfoList().empty())
            this->Recreate(pPfo, nullptr, newPfoList);
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_recreatedClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_recreatedClusterListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_recreatedVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_recreatedVertexListName));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, m_recreatedPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, m_recreatedPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *MasterAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const CaloHit *const pParentCaloHit) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void*>(pParentCaloHit);

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *MasterAlgorithm::CreateCluster(const Cluster *const pInputCluster, const CaloHitList &newCaloHitList,
    const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *MasterAlgorithm::CreateVertex(const Vertex *const pInputVertex) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition();
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *MasterAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ClusterList &newClusterList,
    const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;
    parameters.m_propertiesToAdd = pInputPfo->GetPropertiesMap();

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPC &larTPC, const DetectorGapList &gapList, const std::string &settingsFile, const std::string &name) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora(name));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterCustomContent(pPandora));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The LArTPC
    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = larTPC.GetLArTPCVolumeId();
    larTPCParameters.m_centerX = larTPC.GetCenterX();
    larTPCParameters.m_centerY = larTPC.GetCenterY();
    larTPCParameters.m_centerZ = larTPC.GetCenterZ();
    larTPCParameters.m_widthX = larTPC.GetWidthX();
    larTPCParameters.m_widthY = larTPC.GetWidthY();
    larTPCParameters.m_widthZ = larTPC.GetWidthZ();
    larTPCParameters.m_wirePitchU = larTPC.GetWirePitchU();
    larTPCParameters.m_wirePitchV = larTPC.GetWirePitchV();
    larTPCParameters.m_wirePitchW = larTPC.GetWirePitchW();
    larTPCParameters.m_wireAngleU = larTPC.GetWireAngleU();
    larTPCParameters.m_wireAngleV = larTPC.GetWireAngleV();
    larTPCParameters.m_wireAngleW = larTPC.GetWireAngleW();
    larTPCParameters.m_sigmaUVW = larTPC.GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = larTPC.IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    const float tpcMinX(larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX()), tpcMaxX(larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX());

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap && (((pLineGap->GetLineEndX() >= tpcMinX) && (pLineGap->GetLineEndX() <= tpcMaxX)) ||
            ((pLineGap->GetLineStartX() >= tpcMinX) && (pLineGap->GetLineStartX() <= tpcMaxX))) )
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            const LineGapType lineGapType(pLineGap->GetLineGapType());
            lineGapParameters.m_lineGapType = lineGapType;
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();

            if (m_fullWidthCRWorkerWireGaps && ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W)))
            {
                lineGapParameters.m_lineStartX = -std::numeric_limits<float>::max();
                lineGapParameters.m_lineEndX = std::numeric_limits<float>::max();
            }

            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPCMap &larTPCMap, const DetectorGapList &gapList, const std::string &settingsFile, const std::string &name) const
{
    if (larTPCMap.empty())
    {
        std::cout << "MasterAlgorithm::CreateWorkerInstance - no LArTPC details provided" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    // The Pandora instance
    const Pandora *const pPandora(new Pandora(name));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RegisterCustomContent(pPandora));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The Parent LArTPC
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = 0;
    larTPCParameters.m_centerX = 0.5f * (parentMaxX + parentMinX);
    larTPCParameters.m_centerY = 0.5f * (parentMaxY + parentMinY);
    larTPCParameters.m_centerZ = 0.5f * (parentMaxZ + parentMinZ);
    larTPCParameters.m_widthX = parentMaxX - parentMinX;
    larTPCParameters.m_widthY = parentMaxY - parentMinY;
    larTPCParameters.m_widthZ = parentMaxZ - parentMinZ;
    larTPCParameters.m_wirePitchU = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchV = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchW = pFirstLArTPC->GetWirePitchW();
    larTPCParameters.m_wireAngleU = pFirstLArTPC->GetWireAngleU();
    larTPCParameters.m_wireAngleV = pFirstLArTPC->GetWireAngleV();
    larTPCParameters.m_wireAngleW = pFirstLArTPC->GetWireAngleW();
    larTPCParameters.m_sigmaUVW = pFirstLArTPC->GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = pFirstLArTPC->IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap)
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            lineGapParameters.m_lineGapType = pLineGap->GetLineGapType();
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();
            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RegisterCustomContent(const Pandora *const /*pPandora*/) const
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    ExternalSteeringParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        pExternalParameters = dynamic_cast<ExternalSteeringParameters*>(this->GetExternalParameters());
        if (!pExternalParameters) return STATUS_CODE_FAILURE;
    }

    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "SliceSelectionTools",
            algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            SliceSelectionBaseTool *const pSliceSelectionTool(dynamic_cast<SliceSelectionBaseTool*>(pAlgorithmTool));
            if (!pSliceSelectionTool) return STATUS_CODE_INVALID_PARAMETER;
            m_sliceSelectionToolVector.push_back(pSliceSelectionTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunAllHitsCosmicReco, xmlHandle, "ShouldRunAllHitsCosmicReco", m_shouldRunAllHitsCosmicReco));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunStitching, xmlHandle, "ShouldRunStitching", m_shouldRunStitching));

    if (m_shouldRunStitching && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunStitching requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunStitching)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "StitchingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            StitchingBaseTool *const pStitchingTool(dynamic_cast<StitchingBaseTool*>(pAlgorithmTool));
            if (!pStitchingTool) return STATUS_CODE_INVALID_PARAMETER;
            m_stitchingToolVector.push_back(pStitchingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicHitRemoval, xmlHandle, "ShouldRunCosmicHitRemoval", m_shouldRunCosmicHitRemoval));

    if (m_shouldRunCosmicHitRemoval && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunCosmicHitRemoval requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "CosmicRayTaggingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool(dynamic_cast<CosmicRayTaggingBaseTool*>(pAlgorithmTool));
            if (!pCosmicRayTaggingTool) return STATUS_CODE_INVALID_PARAMETER;
            m_cosmicRayTaggingToolVector.push_back(pCosmicRayTaggingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunSlicing, xmlHandle, "ShouldRunSlicing", m_shouldRunSlicing));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunNeutrinoRecoOption, xmlHandle, "ShouldRunNeutrinoRecoOption", m_shouldRunNeutrinoRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicRecoOption, xmlHandle, "ShouldRunCosmicRecoOption", m_shouldRunCosmicRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldPerformSliceId, xmlHandle, "ShouldPerformSliceId", m_shouldPerformSliceId));

    if (m_shouldPerformSliceId && (!m_shouldRunSlicing || !m_shouldRunNeutrinoRecoOption || !m_shouldRunCosmicRecoOption))
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldPerformSliceId requires ShouldRunSlicing and both neutrino and cosmic reconstruction options" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldPerformSliceId)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "SliceIdTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            SliceIdBaseTool *const pSliceIdIdTool(dynamic_cast<SliceIdBaseTool*>(pAlgorithmTool));
            if (!pSliceIdIdTool) return STATUS_CODE_INVALID_PARAMETER;
            m_sliceIdToolVector.push_back(pSliceIdIdTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_printOverallRecoStatus, xmlHandle, "PrintOverallRecoStatus", m_printOverallRecoStatus));


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeOverallRecoStatus", m_visualizeOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LArCaloHitVersion", m_larCaloHitVersion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRemoveOutOfTimeHits", m_shouldRemoveOutOfTimeHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullWidthCRWorkerWireGaps", m_fullWidthCRWorkerWireGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PassMCParticlesToWorkerInstances", m_passMCParticlesToWorkerInstances));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));
 

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CRSettingsFile", m_crSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuSettingsFile", m_nuSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SlicingSettingsFile", m_slicingSettingsFile));
    m_crSettingsFile = LArFileHelper::FindFileInPath(m_crSettingsFile, m_filePathEnvironmentVariable);
    m_nuSettingsFile = LArFileHelper::FindFileInPath(m_nuSettingsFile, m_filePathEnvironmentVariable);
    m_slicingSettingsFile = LArFileHelper::FindFileInPath(m_slicingSettingsFile, m_filePathEnvironmentVariable);

    if (m_passMCParticlesToWorkerInstances)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputMCParticleListName", m_inputMCParticleListName));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputHitListName", m_inputHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedPfoListName", m_recreatedPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedClusterListName", m_recreatedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedVertexListName", m_recreatedVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMaxX0", m_inTimeMaxX0));



    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PfoListName", m_pfoListName));
 


    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "TrackDirection", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
      {
        TrackDirectionBaseTool *const pTrackDirectionTool(dynamic_cast<TrackDirectionBaseTool*>(pAlgorithmTool));
	if (!pTrackDirectionTool) return STATUS_CODE_INVALID_PARAMETER;
	m_trackDirectionToolVector.push_back(pTrackDirectionTool);
      }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const InputBool inputBool,
    const TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool)
{
    if (pExternalParameters && inputBool.IsInitialized())
    {
        outputBool = inputBool.Get();
    }
    else
    {
  
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, xmlTag, outputBool));
	
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
