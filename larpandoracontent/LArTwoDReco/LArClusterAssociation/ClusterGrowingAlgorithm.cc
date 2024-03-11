/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterGrowingAlgorithm::ClusterGrowingAlgorithm() :
    m_maxClusterSeparation(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterGrowingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;

    if (m_inputClusterListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));
    }

    if (!pClusterList || pClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ClusterGrowingAlgorithm: unable to find cluster list " << m_inputClusterListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    ClusterVector inputClusters, seedClusters;
    this->GetListOfCleanClusters(pClusterList, inputClusters);
    this->GetListOfSeedClusters(inputClusters, seedClusters);

    while (true)
    {
        ClusterVector currentClusters, nonSeedClusters;
        this->GetListOfCleanClusters(pClusterList, currentClusters);
        this->GetListOfNonSeedClusters(currentClusters, seedClusters, nonSeedClusters);

        ClusterMergeMap clusterMergeMap;
        this->PopulateClusterMergeMap(seedClusters, nonSeedClusters, clusterMergeMap);

        if (clusterMergeMap.empty())
            break;

        this->MergeClusters(clusterMergeMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterGrowingAlgorithm::GetListOfNonSeedClusters(
    const ClusterVector &inputClusters, const ClusterVector &seedClusters, ClusterVector &nonSeedClusters) const
{
    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (seedClusters.end() != std::find(seedClusters.begin(), seedClusters.end(), pCluster))
            continue;

        nonSeedClusters.push_back(pCluster);
    }

    std::sort(nonSeedClusters.begin(), nonSeedClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterGrowingAlgorithm::PopulateClusterMergeMap(
    const ClusterVector &seedClusters, const ClusterVector &nonSeedClusters, ClusterMergeMap &clusterMergeMap) const
{
    for (ClusterVector::const_iterator nIter = nonSeedClusters.begin(), nIterEnd = nonSeedClusters.end(); nIter != nIterEnd; ++nIter)
    {
        const Cluster *const pNonSeedCluster = *nIter;

        const Cluster *pBestSeedCluster(NULL);
        float bestDistance(m_maxClusterSeparation);

        for (ClusterVector::const_iterator sIter = seedClusters.begin(), sIterEnd = seedClusters.end(); sIter != sIterEnd; ++sIter)
        {
            const Cluster *const pThisSeedCluster = *sIter;
            const float thisDistance(LArClusterHelper::GetClosestDistance(pNonSeedCluster, pThisSeedCluster));

            if (thisDistance < bestDistance)
            {
                pBestSeedCluster = pThisSeedCluster;
                bestDistance = thisDistance;
            }
        }

        if (pBestSeedCluster)
            clusterMergeMap[pBestSeedCluster].push_back(pNonSeedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterGrowingAlgorithm::MergeClusters(const ClusterMergeMap &clusterMergeMap) const
{
    ClusterList parentClusterList;
    for (const auto &mapEntry : clusterMergeMap)
        parentClusterList.push_back(mapEntry.first);
    parentClusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : parentClusterList)
    {
        const ClusterList &clusterList(clusterMergeMap.at(pParentCluster));

        if (clusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        for (const Cluster *const pAssociatedCluster : clusterList)
        {
            if (m_inputClusterListName.empty())
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pAssociatedCluster));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                    PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pAssociatedCluster, m_inputClusterListName, m_inputClusterListName));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
