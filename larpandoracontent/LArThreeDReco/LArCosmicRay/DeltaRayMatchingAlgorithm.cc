/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray shower matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

DeltaRayMatchingAlgorithm::DeltaRayMatchingAlgorithm() :
    m_minCaloHitsPerCluster(3),
    m_xOverlapWindow(1.f),
    m_distanceForMatching(5.f),
    m_pseudoChi2Cut(3.f),
    m_searchRegion1D(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMatchingAlgorithm::Run()
{
    PfoVector pfoVector;
    this->GetAllPfos(m_parentPfoListName, pfoVector);

    if (pfoVector.empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "DeltaRayMatchingAlgorithm: pfo list " << m_parentPfoListName << " unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    this->InitializeNearbyClusterMaps();

    ClusterLengthMap clusterLengthMap;
    this->ThreeViewMatching(clusterLengthMap);
    this->TwoViewMatching(clusterLengthMap);
    this->OneViewMatching(clusterLengthMap);

    this->ClearNearbyClusterMaps();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::InitializeNearbyClusterMaps()
{
    this->ClearNearbyClusterMaps();
    this->InitializeNearbyClusterMap(m_inputClusterListNameU, m_nearbyClustersU);
    this->InitializeNearbyClusterMap(m_inputClusterListNameV, m_nearbyClustersV);
    this->InitializeNearbyClusterMap(m_inputClusterListNameW, m_nearbyClustersW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::InitializeNearbyClusterMap(const std::string &clusterListName, ClusterToClustersMap &nearbyClusters)
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList))

    if ((NULL == pClusterList) || pClusterList->empty())
        return;

    HitToClusterMap hitToClusterMap;
    CaloHitList allCaloHits;

    for (const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);
        allCaloHits.insert(allCaloHits.end(), daughterHits.begin(), daughterHits.end());

        for (const CaloHit *const pCaloHit : daughterHits)
            (void)hitToClusterMap.insert(HitToClusterMap::value_type(pCaloHit, pCluster));
    }

    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    for (const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);

        for (const CaloHit *const pCaloHit : daughterHits)
        {
            KDTreeBox searchRegionHits = build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D);

            HitKDNode2DList found;
            kdTree.search(searchRegionHits, found);

            for (const auto &hit : found)
            {
                ClusterList &nearbyClusterList(nearbyClusters[pCluster]);
                const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));

                if (nearbyClusterList.end() == std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster))
                    nearbyClusterList.push_back(pNearbyCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::ClearNearbyClusterMaps()
{
    m_nearbyClustersU.clear();
    m_nearbyClustersV.clear();
    m_nearbyClustersW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetAllPfos(const std::string &inputPfoListName, PfoVector &pfoVector) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
        pfoVector.push_back(*iter);

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetTrackPfos(const std::string &inputPfoListName, PfoVector &pfoVector) const
{
    PfoVector inputVector;
    this->GetAllPfos(inputPfoListName, inputVector);

    for (PfoVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pPfo = *iter;

        if (!LArPfoHelper::IsTrack(pPfo))
            continue;

        pfoVector.push_back(pPfo);
    }

    // ATTN Track pfo list is sorted only because the inputVector is sorted
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetClusters(const std::string &inputClusterListName, pandora::ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, inputClusterListName, pClusterList))

    if (NULL == pClusterList)
        return;

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (!pCluster->IsAvailable() || (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::ThreeViewMatching(ClusterLengthMap &clusterLengthMap) const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    PfoLengthMap pfoLengthMap;
    ParticleList initialParticleList, finalParticleList;
    this->ThreeViewMatching(clustersU, clustersV, clustersW, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->SelectParticles(initialParticleList, clusterLengthMap, finalParticleList);
    this->CreateParticles(finalParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::TwoViewMatching(ClusterLengthMap &clusterLengthMap) const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    PfoLengthMap pfoLengthMap;
    ParticleList initialParticleList, finalParticleList;
    this->TwoViewMatching(clustersU, clustersV, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->TwoViewMatching(clustersV, clustersW, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->TwoViewMatching(clustersW, clustersU, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->SelectParticles(initialParticleList, clusterLengthMap, finalParticleList);
    this->CreateParticles(finalParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::OneViewMatching(ClusterLengthMap &clusterLengthMap) const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    PfoLengthMap pfoLengthMap;
    ParticleList initialParticleList, finalParticleList;
    this->ThreeViewMatching(clustersU, clustersV, clustersW, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->OneViewMatching(clustersU, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->OneViewMatching(clustersV, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->OneViewMatching(clustersW, clusterLengthMap, pfoLengthMap, initialParticleList);
    this->SelectParticles(initialParticleList, clusterLengthMap, finalParticleList);
    this->CreateParticles(finalParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::ThreeViewMatching(const ClusterVector &clusters1, const ClusterVector &clusters2,
    const ClusterVector &clusters3, ClusterLengthMap &clusterLengthMap, PfoLengthMap &pfoLengthMap, ParticleList &particleList) const
{
    if (clusters1.empty() || clusters2.empty() || clusters3.empty())
        return;

    for (const Cluster *const pCluster1 : clusters1)
    {
        if (!pCluster1->IsAvailable())
            continue;

        for (const Cluster *const pCluster2 : clusters2)
        {
            if (!pCluster2->IsAvailable())
                continue;

            for (const Cluster *const pCluster3 : clusters3)
            {
                if (!pCluster3->IsAvailable())
                    continue;

                if (!this->AreClustersMatched(pCluster1, pCluster2, pCluster3))
                    continue;

                const ParticleFlowObject *pBestPfo = NULL;
                this->FindBestParentPfo(pCluster1, pCluster2, pCluster3, clusterLengthMap, pfoLengthMap, pBestPfo);

                // ATTN Need to record all matches when all three views are used
                particleList.push_back(Particle(pCluster1, pCluster2, pCluster3, pBestPfo));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::TwoViewMatching(const ClusterVector &clusters1, const ClusterVector &clusters2,
    ClusterLengthMap &clusterLengthMap, PfoLengthMap &pfoLengthMap, ParticleList &particleList) const
{
    if (clusters1.empty() || clusters2.empty())
        return;

    for (const Cluster *const pCluster1 : clusters1)
    {
        if (!pCluster1->IsAvailable())
            continue;

        for (const Cluster *const pCluster2 : clusters2)
        {
            if (!pCluster2->IsAvailable())
                continue;

            if (!this->AreClustersMatched(pCluster1, pCluster2, NULL))
                continue;

            const ParticleFlowObject *pBestPfo = NULL;
            this->FindBestParentPfo(pCluster1, pCluster2, NULL, clusterLengthMap, pfoLengthMap, pBestPfo);

            if (NULL == pBestPfo)
                continue;

            particleList.push_back(Particle(pCluster1, pCluster2, NULL, pBestPfo));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::OneViewMatching(
    const ClusterVector &clusters, ClusterLengthMap &clusterLengthMap, PfoLengthMap &pfoLengthMap, ParticleList &particleList) const
{
    if (clusters.empty())
        return;

    for (const Cluster *const pCluster : clusters)
    {
        if (!pCluster->IsAvailable())
            continue;

        const ParticleFlowObject *pBestPfo = NULL;
        this->FindBestParentPfo(pCluster, NULL, NULL, clusterLengthMap, pfoLengthMap, pBestPfo);

        if (NULL == pBestPfo)
            continue;

        particleList.push_back(Particle(pCluster, NULL, NULL, pBestPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::SelectParticles(const ParticleList &initialParticles, ClusterLengthMap &clusterLengthMap, ParticleList &finalParticles) const
{
    for (const Particle &particle1 : initialParticles)
    {
        bool isGoodParticle(true);

        for (const Particle &particle2 : initialParticles)
        {
            const bool commonU(particle1.GetClusterU() == particle2.GetClusterU());
            const bool commonV(particle1.GetClusterV() == particle2.GetClusterV());
            const bool commonW(particle1.GetClusterW() == particle2.GetClusterW());

            const bool ambiguousU(commonU && NULL != particle1.GetClusterU());
            const bool ambiguousV(commonV && NULL != particle1.GetClusterV());
            const bool ambiguousW(commonW && NULL != particle1.GetClusterW());

            if (commonU && commonV && commonW)
                continue;

            if (ambiguousU || ambiguousV || ambiguousW)
            {
                if (particle2.GetNViews() > particle1.GetNViews())
                {
                    isGoodParticle = false;
                }
                else if (particle2.GetNViews() == particle1.GetNViews() && NULL != particle2.GetParentPfo())
                {
                    if ((NULL == particle1.GetParentPfo()) || (particle2.GetNCaloHits() > particle1.GetNCaloHits()) ||
                        (particle2.GetNCaloHits() == particle1.GetNCaloHits() &&
                            this->GetLength(particle2, clusterLengthMap) >= this->GetLength(particle1, clusterLengthMap)))
                    {
                        isGoodParticle = false;
                    }
                }

                if (!isGoodParticle)
                    break;
            }
        }

        if (isGoodParticle && NULL != particle1.GetParentPfo())
            finalParticles.push_back(particle1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::CreateParticles(const ParticleList &particleList) const
{
    PfoVector parentVector, daughterVector;
    this->GetTrackPfos(m_parentPfoListName, parentVector);
    this->GetAllPfos(m_daughterPfoListName, daughterVector);

    PfoList parentList(parentVector.begin(), parentVector.end());
    PfoList daughterList(daughterVector.begin(), daughterVector.end());

    for (const Particle &particle : particleList)
    {
        const ParticleFlowObject *const pParentPfo = particle.GetParentPfo();

        if (NULL == pParentPfo)
            continue;

        const Cluster *const pClusterU = particle.GetClusterU();
        const Cluster *const pClusterV = particle.GetClusterV();
        const Cluster *const pClusterW = particle.GetClusterW();

        if (NULL == pClusterU && NULL == pClusterV && NULL == pClusterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterList clusterList;

        if (pClusterU)
            clusterList.push_back(pClusterU);

        if (pClusterV)
            clusterList.push_back(pClusterV);

        if (pClusterW)
            clusterList.push_back(pClusterW);

        if (parentList.end() != std::find(parentList.begin(), parentList.end(), pParentPfo))
        {
            this->CreateDaughterPfo(clusterList, pParentPfo);
        }
        else if (daughterList.end() != std::find(daughterList.begin(), daughterList.end(), pParentPfo))
        {
            this->AddToDaughterPfo(clusterList, pParentPfo);
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMatchingAlgorithm::AreClustersMatched(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3) const
{
    if (nullptr == pCluster1 && nullptr == pCluster2 && nullptr == pCluster3)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // First step: Check X overlap
    float xMin1(-std::numeric_limits<float>::max()), xMax1(+std::numeric_limits<float>::max());
    float xMin2(-std::numeric_limits<float>::max()), xMax2(+std::numeric_limits<float>::max());
    float xMin3(-std::numeric_limits<float>::max()), xMax3(+std::numeric_limits<float>::max());

    if (nullptr != pCluster1)
        pCluster1->GetClusterSpanX(xMin1, xMax1);

    if (nullptr != pCluster2)
        pCluster2->GetClusterSpanX(xMin2, xMax2);

    if (nullptr != pCluster3)
        pCluster3->GetClusterSpanX(xMin3, xMax3);

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, std::max(xMin2, xMin3)) - xPitch);
    const float xMax(std::min(xMax1, std::min(xMax2, xMax3)) + xPitch);
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return false;

    if (nullptr == pCluster1 || nullptr == pCluster2 || nullptr == pCluster3)
        return true;

    // Second step: Check 3D matching
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));
    const float pitch1{LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType1)};
    const float pitch2{LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType2)};
    const float pitch3{LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType3)};
    const float pitchMax{std::max({pitch1, pitch2, pitch3})};

    if (hitType1 == hitType2 || hitType2 == hitType3 || hitType3 == hitType1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const unsigned int nSamplingPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    for (unsigned int n = 0; n < nSamplingPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nSamplingPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMin1(0.f), zMin2(0.f), zMin3(0.f), zMax1(0.f), zMax2(0.f), zMax3(0.f);
            pCluster1->GetClusterSpanZ(xmin, xmax, zMin1, zMax1);
            pCluster2->GetClusterSpanZ(xmin, xmax, zMin2, zMax2);
            pCluster3->GetClusterSpanZ(xmin, xmax, zMin3, zMax3);

            const float z1(0.5f * (zMin1 + zMax1));
            const float z2(0.5f * (zMin2 + zMax2));
            const float z3(0.5f * (zMin3 + zMax3));

            const float dz1(zMax1 - zMin1);
            const float dz2(zMax2 - zMin2);
            const float dz3(zMax3 - zMin3);
            const float dz4(pitchMax);

            const float zproj1(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, z2, z3));
            const float zproj2(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType3, hitType1, z3, z1));
            const float zproj3(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, z1, z2));

            const float deltaSquared(((z1 - zproj1) * (z1 - zproj1) + (z2 - zproj2) * (z2 - zproj2) + (z3 - zproj3) * (z3 - zproj3)) / 3.f);
            const float sigmaSquared(dz1 * dz1 + dz2 * dz2 + dz3 * dz3 + dz4 * dz4);
            const float pseudoChi2(deltaSquared / sigmaSquared);

            if (pseudoChi2 < m_pseudoChi2Cut)
                return true;
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::FindBestParentPfo(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const Cluster *const pCluster3, ClusterLengthMap &clusterLengthMap, PfoLengthMap &pfoLengthMap, const ParticleFlowObject *&pBestPfo) const
{
    PfoVector pfoVector;
    this->GetTrackPfos(m_parentPfoListName, pfoVector);
    this->GetAllPfos(m_daughterPfoListName, pfoVector);

    if (pfoVector.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    unsigned int numViews(0);
    float lengthSquared(0.f);

    if (pCluster1)
    {
        lengthSquared += this->GetLengthFromCache(pCluster1, clusterLengthMap);
        ++numViews;
    }

    if (pCluster2)
    {
        lengthSquared += this->GetLengthFromCache(pCluster2, clusterLengthMap);
        ++numViews;
    }

    if (pCluster3)
    {
        lengthSquared += this->GetLengthFromCache(pCluster3, clusterLengthMap);
        ++numViews;
    }

    float bestDistanceSquared(static_cast<float>(numViews) * m_distanceForMatching * m_distanceForMatching);

    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        if (lengthSquared > this->GetLengthFromCache(pPfo, pfoLengthMap))
            continue;

        try
        {
            float distanceSquared(0.f);

            if (NULL != pCluster1)
                distanceSquared += this->GetDistanceSquaredToPfo(pCluster1, pPfo);

            if (NULL != pCluster2)
                distanceSquared += this->GetDistanceSquaredToPfo(pCluster2, pPfo);

            if (NULL != pCluster3)
                distanceSquared += this->GetDistanceSquaredToPfo(pCluster3, pPfo);

            if (distanceSquared < bestDistanceSquared)
            {
                pBestPfo = pPfo;
                bestDistanceSquared = distanceSquared;
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode()))
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::GetLengthFromCache(const Cluster *const pCluster, ClusterLengthMap &clusterLengthMap) const
{
    ClusterLengthMap::const_iterator iter = clusterLengthMap.find(pCluster);

    if (clusterLengthMap.end() != iter)
        return iter->second;

    const float lengthSquared(LArClusterHelper::GetLengthSquared(pCluster));
    (void)clusterLengthMap.insert(ClusterLengthMap::value_type(pCluster, lengthSquared));
    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::GetLengthFromCache(const ParticleFlowObject *const pPfo, PfoLengthMap &pfoLengthMap) const
{
    PfoLengthMap::const_iterator iter = pfoLengthMap.find(pPfo);

    if (pfoLengthMap.end() != iter)
        return iter->second;

    const float lengthSquared(LArPfoHelper::GetTwoDLengthSquared(pPfo));
    (void)pfoLengthMap.insert(PfoLengthMap::value_type(pPfo, lengthSquared));
    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::GetLength(const Particle &particle, ClusterLengthMap &clusterLengthMap) const
{
    float lengthSquared(0.f);

    if (particle.GetClusterU())
        lengthSquared += this->GetLengthFromCache(particle.GetClusterU(), clusterLengthMap);

    if (particle.GetClusterV())
        lengthSquared += this->GetLengthFromCache(particle.GetClusterV(), clusterLengthMap);

    if (particle.GetClusterW())
        lengthSquared += this->GetLengthFromCache(particle.GetClusterW(), clusterLengthMap);

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::GetDistanceSquaredToPfo(const Cluster *const pCluster, const ParticleFlowObject *const pPfo) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    ClusterList comparisonList;
    const ClusterToClustersMap &nearbyClusters((TPC_VIEW_U == hitType) ? m_nearbyClustersU : (TPC_VIEW_V == hitType) ? m_nearbyClustersV : m_nearbyClustersW);

    if (!nearbyClusters.count(pCluster))
        return std::numeric_limits<float>::max();

    ClusterList pfoClusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, pfoClusterList);

    for (const Cluster *const pPfoCluster : pfoClusterList)
    {
        const ClusterList &clusterList(nearbyClusters.at(pCluster));

        if ((clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pPfoCluster)) &&
            (comparisonList.end() == std::find(comparisonList.begin(), comparisonList.end(), pPfoCluster)))
        {
            comparisonList.push_back(pPfoCluster);
        }
    }

    if (comparisonList.empty())
        return std::numeric_limits<float>::max();

    return LArClusterHelper::GetClosestDistance(pCluster, comparisonList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::CreateDaughterPfo(const ClusterList &clusterList, const ParticleFlowObject *const pParentPfo) const
{
    const PfoList *pPfoList = NULL;
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    // TODO Correct these placeholder parameters
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = E_MINUS; // SHOWER
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusterList;

    const ParticleFlowObject *pDaughterPfo(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pDaughterPfo));

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_daughterPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::AddToDaughterPfo(const ClusterList &clusterList, const ParticleFlowObject *const pParentPfo) const
{
    for (const Cluster *const pDaughterCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const std::string clusterListName(
            (TPC_VIEW_U == hitType) ? m_inputClusterListNameU : (TPC_VIEW_V == hitType) ? m_inputClusterListNameV : m_inputClusterListNameW);

        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pParentPfo, hitType, pfoClusters);

        if (pfoClusters.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const Cluster *const pParentCluster = *(pfoClusters.begin());

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster, clusterListName, clusterListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRayMatchingAlgorithm::Particle::Particle(
    const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3, const ParticleFlowObject *const pPfo) :
    m_pClusterU(NULL),
    m_pClusterV(NULL),
    m_pClusterW(NULL),
    m_pParentPfo(NULL)
{
    const HitType hitType1(NULL != pCluster1 ? LArClusterHelper::GetClusterHitType(pCluster1) : HIT_CUSTOM);
    const HitType hitType2(NULL != pCluster2 ? LArClusterHelper::GetClusterHitType(pCluster2) : HIT_CUSTOM);
    const HitType hitType3(NULL != pCluster3 ? LArClusterHelper::GetClusterHitType(pCluster3) : HIT_CUSTOM);

    m_pClusterU = ((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : (TPC_VIEW_U == hitType3) ? pCluster3 : NULL);
    m_pClusterV = ((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : (TPC_VIEW_V == hitType3) ? pCluster3 : NULL);
    m_pClusterW = ((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : (TPC_VIEW_W == hitType3) ? pCluster3 : NULL);
    m_pParentPfo = pPfo;

    if (NULL == m_pClusterU && NULL == m_pClusterV && NULL == m_pClusterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DeltaRayMatchingAlgorithm::Particle::GetNViews() const
{
    unsigned int numViews(0);

    if (NULL != m_pClusterU)
        numViews += 1;

    if (NULL != m_pClusterV)
        numViews += 1;

    if (NULL != m_pClusterW)
        numViews += 1;

    return numViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DeltaRayMatchingAlgorithm::Particle::GetNCaloHits() const
{
    unsigned int numCaloHits(0);

    if (NULL != m_pClusterU)
        numCaloHits += m_pClusterU->GetNCaloHits();

    if (NULL != m_pClusterV)
        numCaloHits += m_pClusterV->GetNCaloHits();

    if (NULL != m_pClusterW)
        numCaloHits += m_pClusterW->GetNCaloHits();

    return numCaloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentPfoListName", m_parentPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterPfoListName", m_daughterPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OverlapWindow", m_xOverlapWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DistanceForMatching", m_distanceForMatching));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PseudoChi2Cut", m_pseudoChi2Cut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_searchRegion1D));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
