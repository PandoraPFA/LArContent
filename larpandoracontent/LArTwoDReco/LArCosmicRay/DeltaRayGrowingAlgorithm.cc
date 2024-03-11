/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeltaRayGrowingAlgorithm::DeltaRayGrowingAlgorithm() :
    m_minCaloHitsPerCluster(2),
    m_minSeedClusterCaloHits(5),
    m_maxSeedClusterLength(10.f),
    m_maxSeedClusterDisplacement(1.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (!pCluster->IsAvailable() || (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    if (inputClusters.empty())
        return;

    const HitType clusterHitType(LArClusterHelper::GetClusterHitType(inputClusters.front()));

    // Get parent and daughter Pfos
    PfoVector parentPfos, daughterPfos;
    this->GetPfos(m_parentPfoListName, parentPfos);
    this->GetPfos(m_daughterPfoListName, daughterPfos);

    ClusterList parentClusters, daughterClusters;

    for (const Pfo *const pParentPfo : parentPfos)
        LArPfoHelper::GetClusters(pParentPfo, clusterHitType, parentClusters);

    for (const Pfo *const pDaughterPfo : daughterPfos)
        LArPfoHelper::GetClusters(pDaughterPfo, clusterHitType, daughterClusters);

    // Select short parent clusters
    for (const Cluster *const pCluster : parentClusters)
    {
        if (LArClusterHelper::GetLengthSquared(pCluster) < m_maxSeedClusterLength * m_maxSeedClusterLength)
            seedClusters.push_back(pCluster);
    }

    // Select all secondary clusters
    for (const Cluster *const pCluster : daughterClusters)
    {
        seedClusters.push_back(pCluster);
    }

    // Select other possible delta rays
    for (const Cluster *const pCluster : inputClusters)
    {
        if (pCluster->GetNCaloHits() < m_minSeedClusterCaloHits)
            continue;

        const float parentDistance(
            parentClusters.empty() ? std::numeric_limits<float>::max() : LArClusterHelper::GetClosestDistance(pCluster, parentClusters));

        if (parentDistance > m_maxSeedClusterDisplacement)
            continue;

        const float daughterDistance(
            daughterClusters.empty() ? std::numeric_limits<float>::max() : LArClusterHelper::GetClosestDistance(pCluster, daughterClusters));

        if (daughterDistance < m_maxSeedClusterDisplacement)
            continue;

        seedClusters.push_back(pCluster);
    }

    std::sort(seedClusters.begin(), seedClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayGrowingAlgorithm::GetPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, inputPfoListName, pPfoList));

    if (!pPfoList)
        return;

    for (const Pfo *const pPfo : *pPfoList)
        pfoVector.push_back(pPfo);

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentPfoListName", m_parentPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterPfoListName", m_daughterPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSeedClusterCaloHits", m_minSeedClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeedClusterLength", m_maxSeedClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxSeedClusterDisplacement", m_maxSeedClusterDisplacement));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
