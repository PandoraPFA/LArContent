/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingContainers.cc
 *
 *  @brief  Implementation of the delta ray matching containers class.
 *
 *  $Log: $
 */

#include "Pandora/PandoraEnumeratedTypes.h"

#include "Objects/Cluster.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingContainers.h"

using namespace pandora;

namespace lar_content
{

DeltaRayMatchingContainers::DeltaRayMatchingContainers() : m_searchRegion1D(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::FillContainers(const PfoList &inputPfoList, const ClusterList &inputClusterList1,
    const ClusterList &inputClusterList2, const ClusterList &inputClusterList3)
{
    this->FillHitToClusterMap(inputClusterList1);
    this->FillHitToClusterMap(inputClusterList2);
    this->FillHitToClusterMap(inputClusterList3);

    this->FillClusterProximityMap(inputClusterList1);
    this->FillClusterProximityMap(inputClusterList2);
    this->FillClusterProximityMap(inputClusterList3);

    this->FillClusterToPfoMaps(inputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::FillHitToClusterMap(const ClusterList &inputClusterList)
{
    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterMap(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::AddToClusterMap(const Cluster *const pCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U)   ? m_hitToClusterMapU
                                     : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV
                                                               : m_hitToClusterMapW);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
        hitToClusterMap[pCaloHit] = pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::FillClusterProximityMap(const ClusterList &inputClusterList)
{
    if (inputClusterList.empty())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(inputClusterList.front()));

    this->BuildKDTree(hitType);

    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterProximityMap(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::BuildKDTree(const HitType hitType)
{
    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U)   ? m_hitToClusterMapU
                                           : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV
                                                                     : m_hitToClusterMapW);
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);

    CaloHitList allCaloHits;

    for (auto &entry : hitToClusterMap)
        allCaloHits.push_back(entry.first);

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));

    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::AddToClusterProximityMap(const Cluster *const pCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U)   ? m_hitToClusterMapU
                                           : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV
                                                                     : m_hitToClusterMapW);
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U)   ? m_clusterProximityMapU
                                             : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV
                                                                       : m_clusterProximityMapW);

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

            ClusterList &nearbyClusterList(clusterProximityMap[pCluster]);

            if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
                nearbyClusterList.push_back(pNearbyCluster);

            ClusterList &invertedNearbyClusterList(clusterProximityMap[pNearbyCluster]);

            if (std::find(invertedNearbyClusterList.begin(), invertedNearbyClusterList.end(), pCluster) == invertedNearbyClusterList.end())
                invertedNearbyClusterList.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::FillClusterToPfoMaps(const PfoList &inputPfoList)
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
        this->AddClustersToPfoMaps(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::AddClustersToPfoMaps(const ParticleFlowObject *const pPfo)
{
    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pPfo, hitType, pfoClusters);

        ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U)   ? m_clusterToPfoMapU
                                         : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV
                                                                   : m_clusterToPfoMapW);

        for (const Cluster *const pCluster : pfoClusters)
        {
            if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
                continue;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::AddClustersToContainers(const ClusterVector &newClusterVector, const PfoVector &pfoVector)
{
    for (const Cluster *const pNewCluster : newClusterVector)
        this->AddToClusterMap(pNewCluster);

    for (unsigned int i = 0; i < newClusterVector.size(); i++)
    {
        const Cluster *const pNewCluster(newClusterVector.at(i));
        const ParticleFlowObject *const pMuonPfo(pfoVector.at(i));

        this->AddToClusterProximityMap(pNewCluster);

        if (pMuonPfo)
            this->AddClustersToPfoMaps(pMuonPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::RemoveClusterFromContainers(const Cluster *const pDeletedCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pDeletedCluster));
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U)   ? m_hitToClusterMapU
                                     : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV
                                                               : m_hitToClusterMapW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U)   ? m_clusterProximityMapU
                                             : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV
                                                                       : m_clusterProximityMapW);
    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U)   ? m_clusterToPfoMapU
                                     : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV
                                                               : m_clusterToPfoMapW);

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

    const DeltaRayMatchingContainers::ClusterToPfoMap::const_iterator clusterToPfoIter(clusterToPfoMap.find(pDeletedCluster));

    if (clusterToPfoIter != clusterToPfoMap.end())
        clusterToPfoMap.erase(clusterToPfoIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingContainers::ClearContainers()
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

} // namespace lar_content
