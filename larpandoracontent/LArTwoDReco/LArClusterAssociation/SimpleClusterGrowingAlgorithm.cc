/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/SimpleClusterGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the simple cluster growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/SimpleClusterGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SimpleClusterGrowingAlgorithm::SimpleClusterGrowingAlgorithm() :
    m_minCaloHitsPerCluster(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    for (ClusterVector::const_iterator cIter = inputClusters.begin(), cIterEnd = inputClusters.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        seedClusters.push_back(pCluster);
    }

    std::sort(seedClusters.begin(), seedClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
