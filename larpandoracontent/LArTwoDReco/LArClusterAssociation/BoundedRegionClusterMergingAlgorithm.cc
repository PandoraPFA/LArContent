/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/BoundedRegionClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/BoundedRegionClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

BoundedRegionClusterMergingAlgorithm::BoundedRegionClusterMergingAlgorithm() :
    m_xRegionMin(-100.f),
    m_xRegionMax(0.f),
    m_zRegionMin(235.f),
    m_zRegionMax(400.f),
    m_maxDistance(10.f),
    m_minClusterHits(2)
{
}

StatusCode BoundedRegionClusterMergingAlgorithm::Run()
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
            std::cout << "BoundedRegionClusterMergingAlgorithm: unable to find cluster list " << m_inputClusterListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    std::cout << "Running Bounded Region Clustering for cluster list " << m_inputClusterListName << std::endl;
    while (true)
    {
        ClusterVector clusterVector;
        ClusterMergeMap clusterMergeMap;
        this->GetListOfBoundedRegionClusters(pClusterList, clusterVector);
        if (clusterVector.size() < 2)
            break;

        this->PopulateClusterMergeMap(clusterVector, clusterMergeMap);
        if (clusterMergeMap.empty())
            break;

        this->MergeClusters(clusterVector, clusterMergeMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedRegionClusterMergingAlgorithm::GetListOfBoundedRegionClusters(
    const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;

        CartesianVector minimumPos(0.f, 0.f, 0.f), maximumPos(0.f, 0.f, 0.f);
        LArClusterHelper::GetClusterBoundingBox(pCluster, minimumPos, maximumPos);
        (void)maximumPos; // To get rid of an unused variable warning

        if (minimumPos.GetX() > m_xRegionMin && minimumPos.GetX() < m_xRegionMax && minimumPos.GetZ() > m_zRegionMin && minimumPos.GetZ() < m_zRegionMax)
            clusterVector.emplace_back(pCluster);
    }

    if (!clusterVector.empty())
        std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedRegionClusterMergingAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector usedClusters;

    for (unsigned int i = 0; i < clusterVector.size() - 1; ++i)
    {
        const Cluster *const pCluster1 = clusterVector.at(i);

        for (unsigned int j = i + 1; j < clusterVector.size(); ++j)
        {
            const Cluster *const pCluster2 = clusterVector.at(j);
            if (!pCluster2->IsAvailable())
                continue;

            if (std::find(usedClusters.begin(), usedClusters.end(), pCluster2) != usedClusters.end())
                continue;

            if (this->AreClustersAssociated(pCluster1, pCluster2))
            {
                clusterMergeMap[pCluster1].push_back(pCluster2);
                usedClusters.emplace_back(pCluster2);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BoundedRegionClusterMergingAlgorithm::AreClustersAssociated(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    const float distance = LArClusterHelper::GetClosestDistance(pCluster1, pCluster2);
    if (distance < m_maxDistance)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BoundedRegionClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumX", m_xRegionMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaximumX", m_xRegionMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumZ", m_zRegionMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaximumZ", m_zRegionMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistance", m_maxDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));
    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
