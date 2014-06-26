/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the cosimic ray shower growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{

void CosmicRayShowerGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable() || (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    if (inputClusters.empty())
        return;

    // Get hit type
    const Cluster* pFirstCluster = *(inputClusters.begin());
    const HitType clusterHitType(LArThreeDHelper::GetClusterHitType(pFirstCluster));

    // Get primary and secondary pfos
    PfoVector primaryPfos, secondaryPfos;
    this->GetPfos(m_primaryPfoListName, primaryPfos);
    this->GetPfos(m_secondaryPfoListName, secondaryPfos);

    // grow showers from primary and secondary Pfos
    if (m_growPfos)
    {
        this->SelectPrimaryPfoSeeds(primaryPfos, clusterHitType, seedClusters);
        this->SelectSecondaryPfoSeeds(secondaryPfos, clusterHitType, seedClusters);
    }

    // grow showers from available clusters
    if (m_growClusters)
        this->SelectClusterSeeds(inputClusters, primaryPfos, clusterHitType, seedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::GetPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
        pfoVector.push_back(*pIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::GetPfoClusters(const PfoVector &pfoVector, const HitType hitType, ClusterList &clusterList) const
{
    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;
        this->GetPfoClusters(pPfo, hitType, clusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::GetPfoClusters(const ParticleFlowObject *pPfo, const HitType hitType, ClusterList &clusterList) const
{
    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pPfoCluster = *cIter;

        if (hitType != LArThreeDHelper::GetClusterHitType(pPfoCluster))
            continue;

        clusterList.insert(pPfoCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::SelectPrimaryPfoSeeds(const PfoVector &primaryPfos, const HitType clusterHitType,
    ClusterVector &seedClusters) const
{
    for (PfoVector::const_iterator pIter = primaryPfos.begin(), pIterEnd = primaryPfos.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;

        if (!pPfo->GetDaughterPfoList().empty())
            continue;

        ClusterList pfoClusterList;
        this->GetPfoClusters(pPfo, clusterHitType, pfoClusterList);

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pPfoCluster = *cIter;

            if (LArClusterHelper::GetLengthSquared(pPfoCluster) > m_maxSeedClusterLength  * m_maxSeedClusterLength)
                continue;

            seedClusters.push_back(pPfoCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::SelectSecondaryPfoSeeds(const PfoVector &secondaryPfos, const HitType clusterHitType,
    ClusterVector &seedClusters) const
{
    ClusterList pfoClusterList;
    this->GetPfoClusters(secondaryPfos, clusterHitType, pfoClusterList);

    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        seedClusters.push_back(*cIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::SelectClusterSeeds(const ClusterVector &inputClusters, const PfoVector &primaryPfos,
    const HitType clusterHitType, ClusterVector &seedClusters) const
{
    ClusterList pfoClusterList;
    this->GetPfoClusters(primaryPfos, clusterHitType, pfoClusterList);

    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (clusterHitType != LArThreeDHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        bool isSeed(false);

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *pPfoCluster = *cIter;

            if (LArClusterHelper::GetClosestDistance(pCluster, pPfoCluster) < m_maxSeedClusterDisplacement)
            {
                isSeed = true;
                break;
            }
        }

        if (isSeed)
            seedClusters.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryPfoListName", m_primaryPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SecondaryPfoListName", m_secondaryPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "GrowPfos", m_growPfos));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "GrowClusters", m_growClusters));

    m_minCaloHitsPerCluster = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_maxSeedClusterLength = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterLength", m_maxSeedClusterLength));

    m_maxSeedClusterDisplacement = 1.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterDisplacement", m_maxSeedClusterDisplacement));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
