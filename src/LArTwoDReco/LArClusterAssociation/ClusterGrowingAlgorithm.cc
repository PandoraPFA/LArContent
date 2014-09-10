/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ClusterGrowingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;

    if (m_inputClusterListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));
    }

    if (NULL == pClusterList)
    {
        std::cout << "ClusterGrowingAlgorithm: could not find cluster list " << m_inputClusterListName << std::endl;
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

// --- BEGIN EVENT DISPLAY ---
// ClusterList tempList1, tempList2;
// tempList1.insert(seedClusters.begin(), seedClusters.end());
// tempList2.insert(nonSeedClusters.begin(), nonSeedClusters.end());
// PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "Seed Clusters", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "NonSeed Clusters", BLUE);
// PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---

        ClusterMergeMap clusterMergeMap;
        this->PopulateClusterMergeMap(seedClusters, nonSeedClusters, clusterMergeMap);

        if (clusterMergeMap.empty())
            break;

        this->MergeClusters(clusterMergeMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterGrowingAlgorithm::GetListOfNonSeedClusters(const ClusterVector &inputClusters, const ClusterVector &seedClusters,
    ClusterVector &nonSeedClusters) const
{
    const ClusterList seedList(seedClusters.begin(), seedClusters.end());

    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (seedList.count(pCluster))
            continue;

        nonSeedClusters.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterGrowingAlgorithm::PopulateClusterMergeMap(const ClusterVector &seedClusters, const ClusterVector &nonSeedClusters,
    ClusterMergeMap &clusterMergeMap) const
{
    for (ClusterVector::const_iterator nIter = nonSeedClusters.begin(), nIterEnd = nonSeedClusters.end(); nIter != nIterEnd; ++nIter)
    {
        Cluster* pNonSeedCluster = *nIter;

        Cluster* pBestSeedCluster(NULL);
        float bestDistance(m_maxClusterSeparation);

        for (ClusterVector::const_iterator sIter = seedClusters.begin(), sIterEnd = seedClusters.end(); sIter != sIterEnd; ++sIter)
        {
            Cluster* pThisSeedCluster = *sIter;
            const float thisDistance(LArClusterHelper::GetClosestDistance(pNonSeedCluster, pThisSeedCluster));

            if (thisDistance < bestDistance)
            {
                pBestSeedCluster = pThisSeedCluster;
                bestDistance = thisDistance;
            }
        }

        if (pBestSeedCluster)
            clusterMergeMap[pBestSeedCluster].insert(pNonSeedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterGrowingAlgorithm::MergeClusters(const ClusterMergeMap &clusterMergeMap) const
{
    for (ClusterMergeMap::const_iterator sIter = clusterMergeMap.begin(), sIterEnd = clusterMergeMap.end(); sIter != sIterEnd; ++sIter)
    {
        const Cluster *pCluster = sIter->first;
        const ClusterList &clusterList = sIter->second;

        if (clusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Cluster *pSeedCluster = const_cast<Cluster*>(pCluster);

        for (ClusterList::const_iterator nIter = clusterList.begin(), nIterEnd = clusterList.end(); nIter != nIterEnd; ++nIter)
        {
            Cluster* pAssociatedCluster = *nIter;

            if (m_inputClusterListName.empty())
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster,
                    m_inputClusterListName, m_inputClusterListName));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_inputClusterListName.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputClusterListName", m_inputClusterListName));

    m_maxClusterSeparation = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterSeparation", m_maxClusterSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
