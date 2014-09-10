/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ClusterMergingAlgorithm::Run()
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
        std::cout << "ClusterMergingAlgorithm: could not find cluster list " << m_inputClusterListName << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    while (true)
    {
        ClusterVector unsortedVector, clusterVector;
        this->GetListOfCleanClusters(pClusterList, unsortedVector);
        this->GetSortedListOfCleanClusters(unsortedVector, clusterVector);

        ClusterMergeMap clusterMergeMap;
        this->PopulateClusterMergeMap(clusterVector, clusterMergeMap);

        if (clusterMergeMap.empty())
            break;

        this->MergeClusters(clusterVector, clusterMergeMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::MergeClusters(ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterVetoMap clusterVetoMap;

    for (ClusterVector::iterator iter1 = clusterVector.begin(), iterEnd1 = clusterVector.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = *iter1;

        ClusterList mergeList;
        this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, clusterMergeMap, clusterVetoMap, mergeList);

        for (ClusterList::iterator iter2 = mergeList.begin(), iterEnd2 = mergeList.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pAssociatedCluster = *iter2;

            if (clusterVetoMap.end() != clusterVetoMap.find(pAssociatedCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if (!pAssociatedCluster->IsAvailable())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if (m_inputClusterListName.empty())
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster,
                    m_inputClusterListName, m_inputClusterListName));
            }

            (void) clusterVetoMap.insert(ClusterVetoMap::value_type(pAssociatedCluster, true));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::CollectAssociatedClusters(Cluster *pSeedCluster, const ClusterMergeMap &clusterMergeMap, ClusterList& associatedClusterList) const
{
    ClusterVetoMap clusterVetoMap;

    this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, clusterMergeMap, clusterVetoMap, associatedClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::CollectAssociatedClusters(Cluster *pSeedCluster, Cluster *pCurrentCluster, const ClusterMergeMap &clusterMergeMap,
    const ClusterVetoMap &clusterVetoMap, ClusterList &associatedClusterList) const
{
    ClusterVetoMap::const_iterator iter0 = clusterVetoMap.find(pCurrentCluster);

    if (iter0 != clusterVetoMap.end())
        return;

    ClusterMergeMap::const_iterator iter1 = clusterMergeMap.find(pCurrentCluster);

    if (iter1 == clusterMergeMap.end())
        return;

    for (ClusterList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        Cluster *pAssociatedCluster = *iter2;

        if (pAssociatedCluster == pSeedCluster)
            continue;

        if (!associatedClusterList.insert(pAssociatedCluster).second)
            continue;

        this->CollectAssociatedClusters(pSeedCluster, pAssociatedCluster, clusterMergeMap, clusterVetoMap, associatedClusterList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::GetSortedListOfCleanClusters(const ClusterVector &inputClusters, ClusterVector &outputClusters) const
{
    ClusterVector pfoClusters, availableClusters;

    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (!pCluster->IsAvailable())
        {
            pfoClusters.push_back(pCluster);
        }
        else
        {
            availableClusters.push_back(pCluster);
        }
    }

    std::sort(pfoClusters.begin(), pfoClusters.end(), LArClusterHelper::SortByNHits);
    std::sort(availableClusters.begin(), availableClusters.end(), LArClusterHelper::SortByNHits);

    outputClusters.insert(outputClusters.end(), pfoClusters.begin(), pfoClusters.end());
    outputClusters.insert(outputClusters.end(), availableClusters.begin(), availableClusters.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_inputClusterListName.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
