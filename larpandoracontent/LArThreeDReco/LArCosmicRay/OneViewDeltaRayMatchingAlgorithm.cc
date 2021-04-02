/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/OneViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the one view delta ray matching algorithm class.
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
  m_overlapExtension(1.f),
  m_minClusterHits(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OneViewDeltaRayMatchingAlgorithm::Run()
{
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        this->FillHitToClusterMap(hitType);
        this->FillClusterProximityMap(hitType);
    }
    
    this->FillClusterToPfoMaps();

    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        this->PerformOneViewMatching(hitType);

    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        this->PerformRecovery(hitType);

    this->ClearContainers();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::FillHitToClusterMap(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));

    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterMap(pCluster);
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

void OneViewDeltaRayMatchingAlgorithm::FillClusterToPfoMaps()
{     
    const PfoList &muonPfoList(this->GetMuonPfoList());

    for (const ParticleFlowObject *const pPfo : muonPfoList)
        this->AddClustersToPfoMaps(pPfo);

    const PfoList &deltaRayPfoList(this->GetDeltaRayPfoList());

    for (const ParticleFlowObject *const pPfo : deltaRayPfoList)
        this->AddClustersToPfoMaps(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------        

void OneViewDeltaRayMatchingAlgorithm::AddClustersToPfoMaps(const ParticleFlowObject *const pPfo)
{
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pPfo, hitType, pfoClusters);

        ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        for (const Cluster *const pCluster : pfoClusters)
        {
            if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
                continue;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }
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
        
        const ClusterProximityMap::const_iterator iter(clusterProximityMap.find(pAvailableCluster));

        // ATTN: Map will update during loop
        if (iter == clusterProximityMap.end())
            continue;
        
        bool found(true);
        const ClusterList &nearbyClusters(clusterProximityMap.at(pAvailableCluster));
        PfoVector nearbyMuonPfoVector;

        while (found)
        {
            found = false;

            const ParticleFlowObject *pClosestMuonPfo(nullptr);
            float closestDistance(std::numeric_limits<float>::max());

            for (const Cluster *const pNearbyCluster : nearbyClusters)
            {
                if (!this->IsMuonPfo(pNearbyCluster))
                    continue;

                if (std::find(nearbyMuonPfoVector.begin(), nearbyMuonPfoVector.end(), clusterToPfoMap.at(pNearbyCluster)) != nearbyMuonPfoVector.end())
                    continue;

                found = true;

                const float separation(LArClusterHelper::GetClosestDistance(pNearbyCluster, pAvailableCluster));
                
                if (separation < closestDistance)
                {
                    closestDistance = separation;
                    pClosestMuonPfo = clusterToPfoMap.at(pNearbyCluster);
                }
            }

            if (pClosestMuonPfo)
                nearbyMuonPfoVector.push_back(pClosestMuonPfo);
        }

        if (nearbyMuonPfoVector.empty())
            continue;
    
        if (this->AddIntoExistingDeltaRay(pAvailableCluster, nearbyMuonPfoVector))
            continue;

        this->CreateDeltaRay(pAvailableCluster, nearbyMuonPfoVector, modifiedClusters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OneViewDeltaRayMatchingAlgorithm::IsMuonPfo(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
    const ClusterToPfoMap::const_iterator iter(clusterToPfoMap.find(pCluster));

    if (iter == clusterToPfoMap.end())
      return false;

    const ParticleFlowObject *const pPfo(iter->second);
    const PfoList &muonPfoList(this->GetMuonPfoList());
    
    return (std::find(muonPfoList.begin(), muonPfoList.end(), pPfo) != muonPfoList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OneViewDeltaRayMatchingAlgorithm::AddIntoExistingDeltaRay(const Cluster *const pAvailableCluster, const PfoVector &nearbyMuonPfoVector)
{
  const HitType &hitType(LArClusterHelper::GetClusterHitType(pAvailableCluster));
  const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
  const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U) ? TPC_VIEW_V : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
  const ClusterToPfoMap &clusterToPfoMap1((projectedHitType1 == TPC_VIEW_U) ? m_clusterToPfoMapU : (projectedHitType1 == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
  const ClusterToPfoMap &clusterToPfoMap2((projectedHitType2 == TPC_VIEW_U) ? m_clusterToPfoMapU : (projectedHitType2 == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

  for (const ParticleFlowObject *const pNearbyMuonPfo : nearbyMuonPfoVector)
  {
      const Cluster *const pProjectedCluster1(this->GetBestProjectedCluster({pAvailableCluster}, pNearbyMuonPfo, projectedHitType1, false));
      const Cluster *const pProjectedCluster2(this->GetBestProjectedCluster({pAvailableCluster}, pNearbyMuonPfo, projectedHitType2, false));

      if ((!pProjectedCluster1) || (!pProjectedCluster2))
          continue;
  
      const ParticleFlowObject *const pPfo1(clusterToPfoMap1.at(pProjectedCluster1));
      const ParticleFlowObject *const pPfo2(clusterToPfoMap2.at(pProjectedCluster2));

      if (pPfo1 == pPfo2)
      {
          PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo1, pAvailableCluster));

          this->AddClustersToPfoMaps(pPfo1);

          return true;
      }
  }

      /*
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
      */

  return false;
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *OneViewDeltaRayMatchingAlgorithm::GetBestProjectedCluster(const ClusterList &deltaRayClusterGroup, const ParticleFlowObject *const pNearbyMuonPfo,
    const HitType &hitType, const bool findAvailable)
{
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pNearbyMuonPfo, hitType, muonClusterList);

    if (muonClusterList.size() != 1)
        return nullptr;

    const ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
    auto muonProximityIter(clusterProximityMap.find(muonClusterList.front()));

    if (muonProximityIter == clusterProximityMap.end())
        return nullptr;

    float spanMinX(0.f), spanMaxX(0.f);
    this->GetClusterSpanX(deltaRayClusterGroup, spanMinX, spanMaxX);
    
    unsigned int highestHit(0);
    const Cluster *pProjectedCluster(nullptr);

    for (const Cluster *const pNearbyCluster : muonProximityIter->second)
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
            pProjectedCluster = pNearbyCluster;
        }
    }

    return pProjectedCluster;
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

void OneViewDeltaRayMatchingAlgorithm::GetClusterSpanX(const ClusterList &clusterList, float &spanMinX, float &spanMaxX)
{
    spanMinX = std::numeric_limits<float>::max();
    spanMaxX = -std::numeric_limits<float>::max();

    for (const Cluster *const pCluster : clusterList)
    {
        float minX(0.f), maxX(0.f);
        pCluster->GetClusterSpanX(minX, maxX);

        if (minX < spanMinX)
            spanMinX = minX;

        if (maxX > spanMaxX)
            spanMaxX = maxX;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::CreateDeltaRay(const Cluster *const pAvailableCluster, const PfoVector &nearbyMuonPfoVector, ClusterSet &modifiedClusters)
{
    ClusterList clusterGroup, consideredClusters;
    this->GetNearbyAvailableClusters(pAvailableCluster, consideredClusters, clusterGroup);

    for (const Cluster *const pModifiedCluster : clusterGroup)
        modifiedClusters.insert(pModifiedCluster);

    const HitType &hitType(LArClusterHelper::GetClusterHitType(pAvailableCluster));
    const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U) ? TPC_VIEW_V : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
    ClusterList projectedClusters1, projectedClusters2;

    for (const ParticleFlowObject *const pNearbyMuonPfo : nearbyMuonPfoVector)
    {
        const Cluster *const pProjectedCluster1(this->GetBestProjectedCluster(clusterGroup, pNearbyMuonPfo, projectedHitType1, true));
        const Cluster *const pProjectedCluster2(this->GetBestProjectedCluster(clusterGroup, pNearbyMuonPfo, projectedHitType2, true));

        if ((!pProjectedCluster1) && (!pProjectedCluster2))
            continue;

        ClusterList consideredClusters1;
        if (pProjectedCluster1)
            this->GetNearbyAvailableClusters(pProjectedCluster1, consideredClusters1, projectedClusters1);

        ClusterList consideredClusters2;
        if (pProjectedCluster2)
            this->GetNearbyAvailableClusters(pProjectedCluster2, consideredClusters2, projectedClusters2);
    }

    const Cluster *const pCluster1(this->MergeClusterGroup(clusterGroup));
    const Cluster *const pCluster2(this->MergeClusterGroup(projectedClusters1));
    const Cluster *const pCluster3(this->MergeClusterGroup(projectedClusters2));

    this->CreatePfos(pCluster1, pCluster2, pCluster3);
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

const Cluster *OneViewDeltaRayMatchingAlgorithm::MergeClusterGroup(const ClusterList &clusterGroup)
{
    if (clusterGroup.empty())
        return nullptr;

    const Cluster *const pClusterToEnlarge(clusterGroup.front());

    if (clusterGroup.size() == 1)
        return pClusterToEnlarge;

    const HitType &hitType(LArClusterHelper::GetClusterHitType(pClusterToEnlarge));
    const std::string inputClusterListName((hitType == TPC_VIEW_U) ? m_inputClusterListNameU : (hitType == TPC_VIEW_V) ? m_inputClusterListNameV : m_inputClusterListNameW);
    
    for (const Cluster *const pClusterToDelete : clusterGroup)
    {
        this->RemoveClusterFromHitContainers(pClusterToDelete);
 
        if (pClusterToDelete != pClusterToEnlarge)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, inputClusterListName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pClusterToDelete));
        }
    }

    this->AddClusterToHitContainers(pClusterToEnlarge);

    return pClusterToEnlarge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::RemoveClusterFromHitContainers(const Cluster *const pDeletedCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pDeletedCluster));
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);

    CaloHitList caloHitList;
    pDeletedCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const HitToClusterMap::const_iterator iter(hitToClusterMap.find(pCaloHit));

        if (iter == hitToClusterMap.end())
            throw StatusCodeException(STATUS_CODE_FAILURE);
        
        hitToClusterMap.erase(iter);
    }

    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pDeletedCluster));

    if (clusterProximityIter != clusterProximityMap.end())
    {
        const ClusterList &nearbyClusterList(clusterProximityIter->second);

        for (const Cluster *const pNearbyCluster : nearbyClusterList)
        {
            const ClusterProximityMap::iterator iter(clusterProximityMap.find(pNearbyCluster));

            if (iter == clusterProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);

            ClusterList::iterator invertedIter(std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pDeletedCluster));
            invertedCloseClusters.erase(invertedIter);
        }
        
        clusterProximityMap.erase(clusterProximityIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::AddClusterToHitContainers(const Cluster *const pNewCluster)
{
    this->AddToClusterMap(pNewCluster);
    this->AddToClusterProximityMap(pNewCluster);
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
        throw StatusCodeException(STATUS_CODE_FAILURE);
 
    const ParticleFlowObject *pPfo(nullptr);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));

    this->AddClustersToPfoMaps(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::PerformRecovery(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    
    for (const Cluster *const pCluster : inputClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;
        
        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;
        
        this->CreatePfos(pCluster, nullptr, nullptr);
    }
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
