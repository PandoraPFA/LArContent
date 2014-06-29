/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the cosimic ray shower growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

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
    const HitType clusterHitType(LArClusterHelper::GetClusterHitType(pFirstCluster));

    // Select seed clusters for growing
    PfoVector inputPfos;
    this->GetPfos(m_inputPfoListName, inputPfos);
    this->SelectSeedClusters(inputClusters, inputPfos, clusterHitType, seedClusters);
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

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::SelectSeedClusters(const ClusterVector &inputClusters, const PfoVector &inputPfos, 
    const HitType &clusterHitType, ClusterVector &seedClusters) const
{
    // Separate primary clusters (main tracks) and secondary clusters (delta rays)
    ClusterVector primaryClusters, secondaryClusters;

    for (PfoVector::const_iterator pIter = inputPfos.begin(), pIterEnd = inputPfos.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;

        if (!pPfo->GetDaughterPfoList().empty())
            continue;

        if (!this->IsPossibleDeltaRay(pPfo, inputPfos))
        {
            LArPfoHelper::GetClusters(pPfo, clusterHitType, primaryClusters);
        }
        else
        {
            LArPfoHelper::GetClusters(pPfo, clusterHitType, secondaryClusters);
        }
    }

    // Select short primary clusters
    for (ClusterVector::const_iterator cIter = primaryClusters.begin(), cIterEnd = primaryClusters.end(); cIter != cIterEnd; ++cIter)
    {  
        Cluster* pCluster = *cIter;

        if (LArClusterHelper::GetLengthSquared(pCluster) > m_maxSeedClusterLength  * m_maxSeedClusterLength)
            continue;

        seedClusters.push_back(pCluster);
    }

    // Select all secondary clusters
    for (ClusterVector::const_iterator cIter = secondaryClusters.begin(), cIterEnd = secondaryClusters.end(); cIter != cIterEnd; ++cIter)
        seedClusters.push_back(*cIter);

    // Select other possible delta rays
    for (ClusterVector::const_iterator cIter = inputClusters.begin(), cIterEnd = inputClusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster* pCluster = *cIter;
 
        if (pCluster->GetNCaloHits() < m_minSeedClusterCaloHits)
            continue;

        const float primaryDistance(LArClusterHelper::GetClosestDistance(pCluster, primaryClusters));
        const float secondaryDistance(LArClusterHelper::GetClosestDistance(pCluster, secondaryClusters));

        if (primaryDistance < m_maxSeedClusterDisplacement && secondaryDistance > m_maxSeedClusterDisplacement)
        {
            seedClusters.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayShowerGrowingAlgorithm::IsPossibleDeltaRay(const ParticleFlowObject *const pDaughterPfo, const PfoVector &pfoVector) const
{
    const float daughterLengthSquared(LArPfoHelper::GetTwoDLengthSquared(pDaughterPfo));

    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pParentPfo = *pIter;

        if (pDaughterPfo == pParentPfo)
            continue;

        if (daughterLengthSquared > 0.5f * LArPfoHelper::GetTwoDLengthSquared(pParentPfo))
            continue;

        if (LArPfoHelper::GetTwoDSeparation(pParentPfo, pDaughterPfo) < m_maxSeedClusterDisplacement)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{   
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    m_minCaloHitsPerCluster = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_minSeedClusterCaloHits = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeedClusterCaloHits", m_minSeedClusterCaloHits));

    m_maxSeedClusterLength = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterLength", m_maxSeedClusterLength));

    m_maxSeedClusterDisplacement = 2.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterDisplacement", m_maxSeedClusterDisplacement));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
