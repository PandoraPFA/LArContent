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

OneViewDeltaRayMatchingAlgorithm::OneViewDeltaRayMatchingAlgorithm() : m_overlapExtension(1.f), m_minClusterHits(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OneViewDeltaRayMatchingAlgorithm::Run()
{
    const PfoList muonPfoList(this->GetMuonPfoList());
    PfoList allPfoList(this->GetDeltaRayPfoList());

    allPfoList.insert(allPfoList.end(), muonPfoList.begin(), muonPfoList.end());

    m_deltaRayMatchingContainers.FillContainers(
        allPfoList, this->GetInputClusterList(TPC_VIEW_U), this->GetInputClusterList(TPC_VIEW_V), this->GetInputClusterList(TPC_VIEW_W));

    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        this->PerformOneViewMatching(hitType);

    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        this->PerformRecovery(hitType);

    m_deltaRayMatchingContainers.ClearContainers();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterList OneViewDeltaRayMatchingAlgorithm::GetInputClusterList(const HitType hitType)
{
    const std::string inputClusterListName((hitType == TPC_VIEW_U)   ? m_inputClusterListNameU
                                           : (hitType == TPC_VIEW_V) ? m_inputClusterListNameV
                                                                     : m_inputClusterListNameW);

    const ClusterList *pInputClusterList(nullptr);

    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, inputClusterListName, pInputClusterList));

    if (!pInputClusterList)
        return ClusterList();

    return *pInputClusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList OneViewDeltaRayMatchingAlgorithm::GetMuonPfoList()
{
    const PfoList *pMuonPfoList(nullptr);

    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_muonPfoListName, pMuonPfoList));

    if (!pMuonPfoList)
        return PfoList();

    return *pMuonPfoList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList OneViewDeltaRayMatchingAlgorithm::GetDeltaRayPfoList()
{
    const PfoList *pDeltaRayPfoList(nullptr);

    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_deltaRayPfoListName, pDeltaRayPfoList));

    if (!pDeltaRayPfoList)
        return PfoList();

    return *pDeltaRayPfoList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::PerformOneViewMatching(const HitType hitType)
{
    const DeltaRayMatchingContainers::ClusterProximityMap &clusterProximityMap(m_deltaRayMatchingContainers.GetClusterProximityMap(hitType));
    const DeltaRayMatchingContainers::ClusterToPfoMap &clusterToPfoMap(m_deltaRayMatchingContainers.GetClusterToPfoMap(hitType));
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

        const DeltaRayMatchingContainers::ClusterProximityMap::const_iterator iter(clusterProximityMap.find(pAvailableCluster));

        // ATTN: Map will update during loop
        if (iter == clusterProximityMap.end())
            continue;

        bool found(false);
        const ClusterList &nearbyClusters(clusterProximityMap.at(pAvailableCluster));
        PfoVector nearbyMuonPfoVector;

        do
        {
            found = false;

            const ParticleFlowObject *pClosestMuonPfo(nullptr);
            float closestDistance(std::numeric_limits<float>::max());

            for (const Cluster *const pNearbyCluster : nearbyClusters)
            {
                if (!this->IsMuonPfo(pNearbyCluster))
                    continue;

                if (std::find(nearbyMuonPfoVector.begin(), nearbyMuonPfoVector.end(), clusterToPfoMap.at(pNearbyCluster)) !=
                    nearbyMuonPfoVector.end())
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
        } while (found);

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
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const DeltaRayMatchingContainers::ClusterToPfoMap &clusterToPfoMap(m_deltaRayMatchingContainers.GetClusterToPfoMap(hitType));
    const DeltaRayMatchingContainers::ClusterToPfoMap::const_iterator iter(clusterToPfoMap.find(pCluster));

    if (iter == clusterToPfoMap.end())
        return false;

    const ParticleFlowObject *const pPfo(iter->second);
    const PfoList &muonPfoList(this->GetMuonPfoList());

    return (std::find(muonPfoList.begin(), muonPfoList.end(), pPfo) != muonPfoList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OneViewDeltaRayMatchingAlgorithm::AddIntoExistingDeltaRay(const Cluster *const pAvailableCluster, const PfoVector &nearbyMuonPfoVector)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pAvailableCluster));
    const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U)   ? TPC_VIEW_V
                                    : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W
                                                                        : TPC_VIEW_U);
    const DeltaRayMatchingContainers::ClusterToPfoMap &clusterToPfoMap1(m_deltaRayMatchingContainers.GetClusterToPfoMap(projectedHitType1));
    const DeltaRayMatchingContainers::ClusterToPfoMap &clusterToPfoMap2(m_deltaRayMatchingContainers.GetClusterToPfoMap(projectedHitType2));

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

            m_deltaRayMatchingContainers.AddClustersToPfoMaps(pPfo1);

            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *OneViewDeltaRayMatchingAlgorithm::GetBestProjectedCluster(
    const ClusterList &deltaRayClusterGroup, const ParticleFlowObject *const pNearbyMuonPfo, const HitType hitType, const bool findAvailable)
{
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pNearbyMuonPfo, hitType, muonClusterList);

    if (muonClusterList.size() != 1)
        return nullptr;

    const DeltaRayMatchingContainers::ClusterProximityMap &clusterProximityMap(m_deltaRayMatchingContainers.GetClusterProximityMap(hitType));
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
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const DeltaRayMatchingContainers::ClusterToPfoMap &clusterToPfoMap(m_deltaRayMatchingContainers.GetClusterToPfoMap(hitType));

    const DeltaRayMatchingContainers::ClusterToPfoMap::const_iterator iter(clusterToPfoMap.find(pCluster));

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

    const HitType hitType(LArClusterHelper::GetClusterHitType(pAvailableCluster));
    const HitType projectedHitType1((hitType == TPC_VIEW_U) ? TPC_VIEW_V : (hitType == TPC_VIEW_V) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType projectedHitType2((projectedHitType1 == TPC_VIEW_U)   ? TPC_VIEW_V
                                    : (projectedHitType1 == TPC_VIEW_V) ? TPC_VIEW_W
                                                                        : TPC_VIEW_U);
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

    this->CreatePfo(pCluster1, pCluster2, pCluster3);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::GetNearbyAvailableClusters(const Cluster *const pCluster, ClusterList &consideredClusters, ClusterList &foundClusters)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const DeltaRayMatchingContainers::ClusterProximityMap &clusterProximityMap(m_deltaRayMatchingContainers.GetClusterProximityMap(hitType));

    consideredClusters.push_back(pCluster);

    if (!pCluster->IsAvailable())
        return;

    if (std::find(foundClusters.begin(), foundClusters.end(), pCluster) == foundClusters.end())
        foundClusters.push_back(pCluster);

    const DeltaRayMatchingContainers::ClusterProximityMap::const_iterator proximityIter(clusterProximityMap.find(pCluster));

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

    const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterToEnlarge));
    const std::string inputClusterListName((hitType == TPC_VIEW_U)   ? m_inputClusterListNameU
                                           : (hitType == TPC_VIEW_V) ? m_inputClusterListNameV
                                                                     : m_inputClusterListNameW);

    for (const Cluster *const pClusterToDelete : clusterGroup)
    {
        m_deltaRayMatchingContainers.RemoveClusterFromContainers(pClusterToDelete);

        if (pClusterToDelete != pClusterToEnlarge)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, inputClusterListName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pClusterToDelete));
        }
    }

    m_deltaRayMatchingContainers.AddClustersToContainers({pClusterToEnlarge}, {nullptr});

    return pClusterToEnlarge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::CreatePfo(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3)
{
    const PfoList *pPfoList(nullptr);
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = E_MINUS;
    pfoParameters.m_charge = PdgTable::GetParticleCharge(E_MINUS);
    pfoParameters.m_mass = PdgTable::GetParticleMass(E_MINUS);
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

    m_deltaRayMatchingContainers.AddClustersToPfoMaps(pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OneViewDeltaRayMatchingAlgorithm::PerformRecovery(const HitType hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));

    for (const Cluster *const pCluster : inputClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;

        this->CreatePfo(pCluster, nullptr, nullptr);
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_deltaRayMatchingContainers.m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OverlapExtension", m_overlapExtension));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
