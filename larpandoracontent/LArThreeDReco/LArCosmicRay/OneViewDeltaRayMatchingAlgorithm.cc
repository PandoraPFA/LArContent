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
    m_availableSearchRegion1D(1.f)
{
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
    const PfoList &muonPfoList(this->GetMuonPfoList());
    
    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);      

    for (const ParticleFlowObject *const pPfo : muonPfoList)
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

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {             
        HitKDNode2DList found;
        const float searchRegion1D(m_availableSearchRegion1D);
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, searchRegion1D, searchRegion1D));

        kdTree.search(searchRegionHits, found);
            
        for (const auto &hit : found)
        {
            const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));
            
            if (pNearbyCluster == pCluster)
                continue;

            const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
            const bool isMuon(clusterToPfoMap.find(pNearbyCluster) != clusterToPfoMap.end());

            if (!isMuon && !pNearbyCluster->IsAvailable())
                continue;
    
            ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? isMuon ? m_muonProximityMapU : m_availableProximityMapU :
                                                     (hitType == TPC_VIEW_V) ? isMuon ? m_muonProximityMapV : m_availableProximityMapV :
                                                     isMuon ? m_muonProximityMapW : m_availableProximityMapW);
  
            ClusterList  &nearbyClusterList(clusterProximityMap[pCluster]);

            if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
                nearbyClusterList.push_back(pNearbyCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::PerformOneViewMatching(const HitType &hitType)
{
    const ClusterProximityMap &muonProximityMap((hitType == TPC_VIEW_U) ? m_muonProximityMapU : (hitType == TPC_VIEW_V) ? m_muonProximityMapV : m_muonProximityMapW);
    const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

    ClusterVector deltaRayClusterList;
    for (auto &entry : muonProximityMap)
    {
        if (entry.first->IsAvailable())
            deltaRayClusterList.push_back(entry.first);
    }

    std::sort(deltaRayClusterList.begin(), deltaRayClusterList.end(), LArClusterHelper::SortByNHits);

    ClusterSet modifiedClusters;
    for (const Cluster *const pDeltaRayCluster : deltaRayClusterList)
    {
        if (modifiedClusters.count(pDeltaRayCluster))
            continue;
        
        // ATTN: Map will update during loop
        const ClusterProximityMap::const_iterator iter(muonProximityMap.find(pDeltaRayCluster));

        if (iter == muonProximityMap.end())
            continue;
        
        const ClusterList &closeMuonClusters(muonProximityMap.at(pDeltaRayCluster));
        const ParticleFlowObject *pClosestMuonPfo(clusterToPfoMap.at(closeMuonClusters.front()));
        
        if (closeMuonClusters.size() != 1)
        {
            float closestDistance(-std::numeric_limits<float>::max());

            for (const Cluster *const pMuonCluster : closeMuonClusters)
            {
                const float separation(LArClusterHelper::GetClosestDistance(pMuonCluster, pDeltaRayCluster));
                
                if (separation < closestDistance)
                {
                    closestDistance = separation;
                    pClosestMuonPfo = clusterToPfoMap.at(pMuonCluster);
                }
            }
        }
        
        // collect close cluster group (send in modified clusters here...)
        ClusterList consideredClusters, clusterGroup;
        this->GetNearbyClusters(pDeltaRayCluster, consideredClusters, clusterGroup);

        for (const Cluster *const pModifiedCluster : clusterGroup)
            modifiedClusters.insert(pModifiedCluster);


        ClusterList projectedClusters1;        
        const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
        this->GetProjectedNearbyClusters(clusterGroup, pClosestMuonPfo, projectedHitType1, projectedClusters1);


        ClusterList projectedClusters2;        
        const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U) ? TPC_VIEW_V : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
        this->GetProjectedNearbyClusters(clusterGroup, pClosestMuonPfo, projectedHitType2, projectedClusters2);

        const Cluster *const pCluster1(this->MergeClusterGroup(clusterGroup));
        const Cluster *const pCluster2(this->MergeClusterGroup(projectedClusters1));
        const Cluster *const pCluster3(this->MergeClusterGroup(projectedClusters2));

        this->CreatePfos(pCluster1, pCluster2, pCluster3);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::GetNearbyClusters(const Cluster *const pCluster, ClusterList &consideredClusters, ClusterList &foundClusters)  
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const ClusterProximityMap &availableProximityMap((hitType == TPC_VIEW_U) ? m_availableProximityMapU : (hitType == TPC_VIEW_V) ? m_availableProximityMapV : m_availableProximityMapW);

    consideredClusters.push_back(pCluster);

    if (!pCluster->IsAvailable())
    {
        std::cout << "ISOBEL THIS SHOULD NOT HAPPEN (1)" << std::endl;
        throw;
    }

    if (std::find(foundClusters.begin(), foundClusters.end(), pCluster) == foundClusters.end())
        foundClusters.push_back(pCluster);
    
    const ClusterProximityMap::const_iterator availableProximityIter(availableProximityMap.find(pCluster));
    
    if (availableProximityIter == availableProximityMap.end())
        return;
    
    for (const Cluster *const pNearbyCluster : availableProximityIter->second)
    {
        if (!pNearbyCluster->IsAvailable())
        {
            std::cout << "ISOBEL THIS SHOULD NOT HAPPEN (2)" << std::endl;
            throw;
        }
        
        if (std::find(consideredClusters.begin(), consideredClusters.end(), pNearbyCluster) != consideredClusters.end())
            continue;

        if (std::find(foundClusters.begin(), foundClusters.end(), pNearbyCluster) != foundClusters.end())
            continue;

        this->GetNearbyClusters(pNearbyCluster, consideredClusters, foundClusters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::GetProjectedNearbyClusters(const ClusterList &deltaRayClusterGroup, const ParticleFlowObject *const pMuonPfo,
    const HitType &hitType, ClusterList &foundClusters)
{
    const ClusterProximityMap &availableProximityMap((hitType == TPC_VIEW_U) ? m_availableProximityMapU : (hitType == TPC_VIEW_V) ? m_availableProximityMapV : m_availableProximityMapW);
    
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pMuonPfo, hitType, muonClusterList);

    if (muonClusterList.size() != 1)
        return;

    if (availableProximityMap.find(muonClusterList.front()) == availableProximityMap.end())
        return;

    float spanMinX(0.f), spanMaxX(0.f);
    this->GetClusterSpanX(deltaRayClusterGroup, spanMinX, spanMaxX);

    const Cluster *pSeedCluster(nullptr); unsigned int highestHit(-std::numeric_limits<unsigned int>::max());
    for (const Cluster *const pNearbyCluster : availableProximityMap.at(muonClusterList.front()))
    {
        float minX(0.f), maxX(0.f);
        pNearbyCluster->GetClusterSpanX(minX, maxX);

        if ((maxX < spanMinX) || (minX > spanMaxX))
            continue;
        
        if (pNearbyCluster->GetNCaloHits() > highestHit)
        {
            highestHit = pNearbyCluster->GetNCaloHits();
            pSeedCluster = pNearbyCluster;
        }
    }

    if (!pSeedCluster)
        return;

    ClusterList consideredClusters;
    this->GetNearbyClusters(pSeedCluster, consideredClusters, foundClusters);
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

    // Remove from the available map
    ClusterProximityMap &availableProximityMap((hitType == TPC_VIEW_U) ? m_availableProximityMapU : (hitType == TPC_VIEW_V) ? m_availableProximityMapV : m_availableProximityMapW);
    const ClusterProximityMap::const_iterator availableProximityIter(availableProximityMap.find(pCluster));

    if (availableProximityIter != availableProximityMap.end())
    {
        const ClusterList &nearbyClusterList(availableProximityIter->second);

        for (const Cluster *const pNearbyCluster : nearbyClusterList)
        {
            const ClusterProximityMap::iterator iter(availableProximityMap.find(pNearbyCluster));

            // ATTN: May have already been deleted (TEST REMOVAL)
            if (iter == availableProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);

            ClusterList::iterator invertedIter(std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pCluster));
            
            if (invertedIter != invertedCloseClusters.end())
                invertedCloseClusters.erase(invertedIter);
        }
        
        availableProximityMap.erase(availableProximityIter);
    }

    // Remove from the muon map
    ClusterProximityMap &muonProximityMap((hitType == TPC_VIEW_U) ? m_muonProximityMapU : (hitType == TPC_VIEW_V) ? m_muonProximityMapV : m_muonProximityMapW);
    const ClusterProximityMap::const_iterator muonProximityIter(muonProximityMap.find(pCluster));

    // shouldn't actually need to do this...
    if (muonProximityIter != muonProximityMap.end())
    {
        const ClusterList &nearbyClusterList(muonProximityIter->second);

        for (const Cluster *const pNearbyCluster : nearbyClusterList)
        {
            const ClusterProximityMap::iterator iter(availableProximityMap.find(pNearbyCluster));

            // ATTN: May have already been deleted (TEST REMOVAL)
            if (iter == muonProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);

            ClusterList::iterator invertedIter(std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pCluster));

            if (invertedIter != invertedCloseClusters.end())
                invertedCloseClusters.erase(invertedIter);
        }
        
        muonProximityMap.erase(muonProximityIter);
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
