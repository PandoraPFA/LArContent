/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/OneViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the unattached delta rays algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/OneViewDeltaRayMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

OneViewDeltaRayMatchingAlgorithm::OneViewDeltaRayMatchingAlgorithm() :
  m_searchRegion1D(2.f),
  m_overlapExtension(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::ClearContainers()
{
    m_hitToClusterMapU.clear();
    m_hitToClusterMapV.clear();
    m_hitToClusterMapW.clear();    

    m_kdTreeU.clear();
    m_kdTreeV.clear();
    m_kdTreeW.clear();    

    m_clusterProximityMapU.clear();
    m_clusterProximityMapV.clear();
    m_clusterProximityMapW.clear();

    m_clusterToPfoMapU.clear();
    m_clusterToPfoMapV.clear();
    m_clusterToPfoMapW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode OneViewDeltaRayMatchingAlgorithm::Run()
{
    this->FillHitToClusterMap(TPC_VIEW_U);
    this->FillHitToClusterMap(TPC_VIEW_V);
    this->FillHitToClusterMap(TPC_VIEW_W);
    
    this->FillClusterToPfoMap(TPC_VIEW_U);
    this->FillClusterToPfoMap(TPC_VIEW_V);
    this->FillClusterToPfoMap(TPC_VIEW_W);
    
    this->FillClusterProximityMap(TPC_VIEW_U);
    this->FillClusterProximityMap(TPC_VIEW_V);
    this->FillClusterProximityMap(TPC_VIEW_W);

    this->PerformOneViewMatching(TPC_VIEW_U);
    this->PerformOneViewMatching(TPC_VIEW_V);
    this->PerformOneViewMatching(TPC_VIEW_W);

    this->PerformRecovery(TPC_VIEW_U);
    this->PerformRecovery(TPC_VIEW_V);
    this->PerformRecovery(TPC_VIEW_W);

    this->ClearContainers();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterList OneViewDeltaRayMatchingAlgorithm::GetInputClusterList(const HitType &hitType)
{
    const std::string inputClusterListName((hitType == TPC_VIEW_U) ? m_inputClusterListNameU : (hitType == TPC_VIEW_V) ? m_inputClusterListNameV : m_inputClusterListNameW);
    
    const ClusterList *pInputClusterList(nullptr);
    
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputClusterListName, pInputClusterList));

    if ((!pInputClusterList) || pInputClusterList->empty())
        return ClusterList();

    return *pInputClusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList OneViewDeltaRayMatchingAlgorithm::GetMuonPfoList()
{
    const PfoList *pMuonPfoList(nullptr);
    
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_muonPfoListName, pMuonPfoList));

    if ((!pMuonPfoList) || pMuonPfoList->empty())
        return PfoList();

    return *pMuonPfoList;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

bool OneViewDeltaRayMatchingAlgorithm::IsMuonPfo(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

    const ClusterToPfoMap::const_iterator iter(clusterToPfoMap.find(pCluster));

    if (iter == clusterToPfoMap.end())
      return false;

    const ParticleFlowObject *const pMuonPfo(iter->second);
    const PfoList &muonPfoList(this->GetMuonPfoList());
    
    return (std::find(muonPfoList.begin(), muonPfoList.end(), pMuonPfo) != muonPfoList.end());
}    
    
//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList OneViewDeltaRayMatchingAlgorithm::GetDeltaRayPfoList()
{
    const PfoList *pDeltaRayPfoList(nullptr);
    
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_deltaRayPfoListName, pDeltaRayPfoList));

    if ((!pDeltaRayPfoList) || pDeltaRayPfoList->empty())
        return PfoList();

    return *pDeltaRayPfoList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OneViewDeltaRayMatchingAlgorithm::IsDeltaRayPfo(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

    const ClusterToPfoMap::const_iterator iter(clusterToPfoMap.find(pCluster));

    if (iter == clusterToPfoMap.end())
      return false;

    const ParticleFlowObject *const pDeltaRayPfo(iter->second);
    const PfoList &deltaRayPfoList(this->GetDeltaRayPfoList());
    
    return (std::find(deltaRayPfoList.begin(), deltaRayPfoList.end(), pDeltaRayPfo) != deltaRayPfoList.end());
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::FillHitToClusterMap(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));

    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterMap(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::AddToClusterMap(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));      
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
        hitToClusterMap[pCaloHit] = pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------        

void OneViewDeltaRayMatchingAlgorithm::FillClusterToPfoMap(const HitType &hitType)
{    
    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);      

    const PfoList &muonPfoList(this->GetMuonPfoList());
    for (const ParticleFlowObject *const pPfo : muonPfoList)
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pPfo, hitType, pfoClusters);
        for (const Cluster *const pCluster : pfoClusters)
            clusterToPfoMap[pCluster] = pPfo;
    }

    const PfoList &deltaRayPfoList(this->GetDeltaRayPfoList());
    for (const ParticleFlowObject *const pPfo : deltaRayPfoList)
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pPfo, hitType, pfoClusters);
        for (const Cluster *const pCluster : pfoClusters)
            clusterToPfoMap[pCluster] = pPfo;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::FillClusterProximityMap(const HitType &hitType)
{
    this->BuildKDTree(hitType);
    
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    
    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterProximityMap(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::BuildKDTree(const HitType &hitType)
{
    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);    

    CaloHitList allCaloHits;
    for (auto &entry : hitToClusterMap)
        allCaloHits.push_back(entry.first);
    
    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));

    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::AddToClusterProximityMap(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));

    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
  
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {             
        HitKDNode2DList found;
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));

        kdTree.search(searchRegionHits, found);
            
        for (const auto &hit : found)
        {
            const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));
            
            if (pNearbyCluster == pCluster)
                continue;
                
            if (std::find(caloHitList.begin(), caloHitList.end(), hit.data) != caloHitList.end())
                continue;
                
            ClusterList  &nearbyClusterList(clusterProximityMap[pCluster]);

            if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
                nearbyClusterList.push_back(pNearbyCluster);
            
            ClusterList &invertedNearbyClusterList(clusterProximityMap[pNearbyCluster]);
            
            if (std::find(invertedNearbyClusterList.begin(), invertedNearbyClusterList.end(), pCluster) == invertedNearbyClusterList.end())
                invertedNearbyClusterList.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::PerformOneViewMatching(const HitType &hitType)
{
    const ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
    const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

    ClusterVector availableClusterList;
    for (auto &entry : clusterProximityMap)
    {
        if (entry.first->IsAvailable())
            availableClusterList.push_back(entry.first);
    }

    std::sort(availableClusterList.begin(), availableClusterList.end(), LArClusterHelper::SortByNHits);

    ClusterSet modifiedClusters;
    for (const Cluster *const pAvailableCluster : availableClusterList)
    {
        if (modifiedClusters.count(pAvailableCluster))
            continue;
        
        // ATTN: Map will update during loop
        const ClusterProximityMap::const_iterator iter(clusterProximityMap.find(pAvailableCluster));

        if (iter == clusterProximityMap.end())
            continue;
        
        const ClusterList &nearbyClusters(clusterProximityMap.at(pAvailableCluster));

        const ParticleFlowObject *pClosestMuonPfo(nullptr);        
	PfoList nearbyMuonPfoList;
	float closestDistance(std::numeric_limits<float>::max());
	for (const Cluster *const pNearbyCluster : nearbyClusters)
        {
	  if (!this->IsMuonPfo(pNearbyCluster))
	    continue;

	  nearbyMuonPfoList.push_back(clusterToPfoMap.at(pNearbyCluster));

	  const float separation(LArClusterHelper::GetClosestDistance(pNearbyCluster, pAvailableCluster));
                
	  if (separation < closestDistance)
          {
	    closestDistance = separation;
	    pClosestMuonPfo = clusterToPfoMap.at(pNearbyCluster);
	  }
        }

	if (!pClosestMuonPfo)
	  continue;
    
        if (this->AddIntoExistingDeltaRay(pAvailableCluster, nearbyMuonPfoList))
	  continue;

	this->CreateDeltaRay(pAvailableCluster, nearbyMuonPfoList, modifiedClusters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::CreateDeltaRay(const Cluster *const pAvailableCluster, const PfoList &nearbyMuonPfoList, ClusterSet &modifiedClusters)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pAvailableCluster));

    // collect close cluster group (send in modified clusters here...)
    ClusterList consideredClusters, clusterGroup;
    this->GetNearbyAvailableClusters(pAvailableCluster, consideredClusters, clusterGroup);

    for (const Cluster *const pModifiedCluster : clusterGroup)
      modifiedClusters.insert(pModifiedCluster);

    ClusterList projectedClusters1;
    const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);

    for (const ParticleFlowObject *const pNearbyMuonPfo : nearbyMuonPfoList)
    {
      this->GetProjectedNearbyClusters(clusterGroup, pNearbyMuonPfo, projectedHitType1, true, projectedClusters1);
    }

    ClusterList projectedClusters2;
    const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U) ? TPC_VIEW_V : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);

    for (const ParticleFlowObject *const pNearbyMuonPfo : nearbyMuonPfoList)
    {
      this->GetProjectedNearbyClusters(clusterGroup, pNearbyMuonPfo, projectedHitType2, true, projectedClusters2);
    }

    const Cluster *const pCluster1(this->MergeClusterGroup(clusterGroup));
    const Cluster *const pCluster2(this->MergeClusterGroup(projectedClusters1));
    const Cluster *const pCluster3(this->MergeClusterGroup(projectedClusters2));

    this->CreatePfos(pCluster1, pCluster2, pCluster3);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OneViewDeltaRayMatchingAlgorithm::AddIntoExistingDeltaRay(const Cluster *const pAvailableCluster, const PfoList &nearbyMuonPfoList)
{
  const HitType &hitType(LArClusterHelper::GetClusterHitType(pAvailableCluster));

  ClusterList projectedClusters1;
  const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
  const ClusterToPfoMap &clusterToPfoMap1((projectedHitType1 == TPC_VIEW_U) ? m_clusterToPfoMapU : (projectedHitType1 == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

  for (const ParticleFlowObject *const pNearbyMuonPfo : nearbyMuonPfoList)
  {
    this->GetProjectedNearbyClusters({pAvailableCluster}, pNearbyMuonPfo, projectedHitType1, false, projectedClusters1);
  }

  ClusterList projectedClusters2;
  const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U) ? TPC_VIEW_V : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
  const ClusterToPfoMap &clusterToPfoMap2((projectedHitType2 == TPC_VIEW_U) ? m_clusterToPfoMapU : (projectedHitType2 == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
  for (const ParticleFlowObject *const pNearbyMuonPfo : nearbyMuonPfoList)
  {
    this->GetProjectedNearbyClusters({pAvailableCluster}, pNearbyMuonPfo, projectedHitType2, false, projectedClusters2);
  }

  for (const Cluster *const pCluster1 : projectedClusters1)
  {
    const ParticleFlowObject *pPfo1(nullptr);
    pPfo1 = clusterToPfoMap1.at(pCluster1);

    for (const Cluster *const pCluster2 : projectedClusters2)
    {
      const ParticleFlowObject *pPfo2(nullptr);
      pPfo2 = clusterToPfoMap2.at(pCluster2);

      if (pPfo1 == pPfo2)
      {
	this->RemoveClusterFromProximityMaps(pAvailableCluster);
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo1, pAvailableCluster));

	return true;
      }
    }
  }

  if ((projectedClusters1.size() == 1) && projectedClusters2.empty())
  {
      const ParticleFlowObject *pPfo1(clusterToPfoMap1.at(projectedClusters1.front()));
      this->RemoveClusterFromProximityMaps(pAvailableCluster);
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo1, pAvailableCluster));

      return true;
  }

  if ((projectedClusters2.size() == 1) && projectedClusters1.empty())
  {
      const ParticleFlowObject *pPfo2(clusterToPfoMap2.at(projectedClusters2.front()));
      this->RemoveClusterFromProximityMaps(pAvailableCluster);
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo2, pAvailableCluster));

      return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::GetNearbyAvailableClusters(const Cluster *const pCluster, ClusterList &consideredClusters, ClusterList &foundClusters)  
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);

    consideredClusters.push_back(pCluster);

    if (!pCluster->IsAvailable())
      return;

    if (std::find(foundClusters.begin(), foundClusters.end(), pCluster) == foundClusters.end())
        foundClusters.push_back(pCluster);
    
    const ClusterProximityMap::const_iterator proximityIter(clusterProximityMap.find(pCluster));
    
    if (proximityIter == clusterProximityMap.end())
        return;
    
    for (const Cluster *const pNearbyCluster : proximityIter->second)
    {
        if (!pNearbyCluster->IsAvailable())
	    continue;
        
        if (std::find(consideredClusters.begin(), consideredClusters.end(), pNearbyCluster) != consideredClusters.end())
            continue;

        if (std::find(foundClusters.begin(), foundClusters.end(), pNearbyCluster) != foundClusters.end())
            continue;

        this->GetNearbyAvailableClusters(pNearbyCluster, consideredClusters, foundClusters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::GetProjectedNearbyClusters(const ClusterList &deltaRayClusterGroup, const ParticleFlowObject *const pMuonPfo,
    const HitType &hitType, const bool findAvailable, ClusterList &foundClusters)
{
    const ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
    
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pMuonPfo, hitType, muonClusterList);

    if (muonClusterList.size() != 1)
        return;

    if (clusterProximityMap.find(muonClusterList.front()) == clusterProximityMap.end())
        return;

    float spanMinX(0.f), spanMaxX(0.f);
    this->GetClusterSpanX(deltaRayClusterGroup, spanMinX, spanMaxX);

    const Cluster *pSeedCluster(nullptr); unsigned int highestHit(-std::numeric_limits<unsigned int>::max());
    for (const Cluster *const pNearbyCluster : clusterProximityMap.at(muonClusterList.front()))
    {
        if (findAvailable && !pNearbyCluster->IsAvailable())
            continue;

        if (!findAvailable && !this->IsDeltaRayPfo(pNearbyCluster))
            continue;

	float minX(0.f), maxX(0.f);
        pNearbyCluster->GetClusterSpanX(minX, maxX);

        if ((maxX < (spanMinX - m_overlapExtension)) || (minX > (spanMaxX + m_overlapExtension)))
            continue;
        
        if (pNearbyCluster->GetNCaloHits() > highestHit)
        {
            highestHit = pNearbyCluster->GetNCaloHits();
            pSeedCluster = pNearbyCluster;
        }
    }

    if (!pSeedCluster)
        return;

    foundClusters.push_back(pSeedCluster);

    // ATTN: Will instantly drop out in the case of findAvailable == false (this is intended)
    ClusterList consideredClusters;
    this->GetNearbyAvailableClusters(pSeedCluster, consideredClusters, foundClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::GetClusterSpanX(const ClusterList &clusterList, float &xMin, float &xMax)
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();

    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const float xCoordinate(pCaloHit->GetPositionVector().GetX());

            if (xCoordinate < xMin)
                xMin = xCoordinate;

            if (xCoordinate > xMax)
                xMax = xCoordinate;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::CreatePfos(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3)
{
    const PfoList *pPfoList(nullptr); std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = E_MINUS;
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);

    if (pCluster1)
        pfoParameters.m_clusterList.push_back(pCluster1);

    if (pCluster2)
        pfoParameters.m_clusterList.push_back(pCluster2);

    if (pCluster3)
        pfoParameters.m_clusterList.push_back(pCluster3);

    if (pfoParameters.m_clusterList.empty())
    {
        std::cout << "YO! ISOBEL THE LIST IS EMPTY" << std::endl;
        throw;
    }

    const ParticleFlowObject *pPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::RemoveClusterFromProximityMaps(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));

    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pCluster));

    if (clusterProximityIter != clusterProximityMap.end())
    {
        const ClusterList &nearbyClusterList(clusterProximityIter->second);

        for (const Cluster *const pNearbyCluster : nearbyClusterList)
        {
            const ClusterProximityMap::iterator iter(clusterProximityMap.find(pNearbyCluster));

            // ATTN: May have already been deleted (TEST REMOVAL)
            if (iter == clusterProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);

            ClusterList::iterator invertedIter(std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pCluster));
	    invertedCloseClusters.erase(invertedIter);
        }
        
        clusterProximityMap.erase(clusterProximityIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *OneViewDeltaRayMatchingAlgorithm::MergeClusterGroup(const ClusterList &clusterGroup)
{
    if (clusterGroup.empty())
        return nullptr;

    const Cluster *const pClusterToEnlarge(clusterGroup.front());
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pClusterToEnlarge));
    const std::string inputClusterListName((hitType == TPC_VIEW_U) ? m_inputClusterListNameU : (hitType == TPC_VIEW_V) ? m_inputClusterListNameV : m_inputClusterListNameW);
    
    for (const Cluster *const pClusterToDelete : clusterGroup)
    {
        this->RemoveClusterFromProximityMaps(pClusterToDelete);
 
        if (pClusterToDelete != pClusterToEnlarge)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, inputClusterListName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pClusterToDelete));
        }
    }

    return pClusterToEnlarge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::PerformRecovery(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    
    for (const Cluster *const pCluster : inputClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;
        
        if (pCluster->GetNCaloHits() < 3)
            continue;
        
        this->CreatePfos(pCluster, nullptr, nullptr);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OneViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName));    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DeltaRayPfoListName", m_deltaRayPfoListName));    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
