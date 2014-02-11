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

namespace lar
{

StatusCode ClusterMergingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    bool carryOn(true);

    while (carryOn)
    {
        carryOn = false;

        ClusterVector clusterVector;
        this->GetListOfCleanClusters(pClusterList, clusterVector);
        std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

        ClusterMergeMap clusterMergeMap;
        this->PopulateClusterMergeMap(clusterVector, clusterMergeMap);

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
                    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));
                (void) clusterVetoMap.insert(ClusterVetoMap::value_type(pAssociatedCluster, true));
                carryOn = true;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
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

StatusCode ClusterMergingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
