/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{

void DeltaRayGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
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

void DeltaRayGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    if (inputClusters.empty())
        return;

    // Get hit type
    const Cluster* pFirstCluster = *(inputClusters.begin());
    const HitType clusterHitType(LArClusterHelper::GetClusterHitType(pFirstCluster));

    // Select seed clusters for growing
    PfoVector inputPfos;
    this->GetPfos(m_inputPfoListName, inputPfos);
    if (inputPfos.empty())
        return;

    this->SelectSeedClusters(inputClusters, inputPfos, clusterHitType, seedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayGrowingAlgorithm::GetPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
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

void DeltaRayGrowingAlgorithm::SelectSeedClusters(const ClusterVector &inputClusters, const PfoVector &inputPfos,
    const HitType clusterHitType, ClusterVector &seedClusters) const
{
    for (ClusterVector::const_iterator iter1 = inputClusters.begin(), iterEnd1 = inputClusters.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pCluster = *iter1;

        if (clusterHitType != LArClusterHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (pCluster->GetNCaloHits() < m_minSeedClusterCaloHits)
            continue;

        bool isSeed(false);

        for (PfoVector::const_iterator iter2 = inputPfos.begin(), iterEnd2 = inputPfos.end(); iter2 != iterEnd2; ++iter2)
        {
            const ParticleFlowObject *pPfo = *iter2;

            ClusterVector pfoClusters;
            LArPfoHelper::GetClusters(pPfo, clusterHitType, pfoClusters);

            if (pfoClusters.empty())
                continue;

            if (LArClusterHelper::GetClosestDistance(pCluster, pfoClusters) < m_maxSeedClusterDisplacement)
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

StatusCode DeltaRayGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{       
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    m_minCaloHitsPerCluster = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_minSeedClusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeedClusterCaloHits", m_minSeedClusterCaloHits));

    m_maxSeedClusterDisplacement = 1.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterDisplacement", m_maxSeedClusterDisplacement));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
