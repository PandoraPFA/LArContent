/**
 *  @file   LArContent/src/LArTwoDSeed/SeedFindingBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the seed finding base algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDSeed/SeedFindingBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedFindingBaseAlgorithm::Run()
{
    // Organize input lists
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    const ClusterList *pNonSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    const bool nonSeedListExists(NULL != pNonSeedClusterList);

    if (!nonSeedListExists)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pNonSeedClusterList));

        if (pSeedClusterList == pNonSeedClusterList)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    // Select the new particle seeds
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pNonSeedClusterList, clusterVector);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    ClusterList seedClusterList;
    this->GetSeedClusterList(clusterVector, seedClusterList);

    // Organize output lists
    if (!seedClusterList.empty())
    {
        if (nonSeedListExists)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListName, m_seedClusterListName, seedClusterList));
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, seedClusterList));
        }
    }

    if (!nonSeedListExists)
    {
        ClusterList nonSeedClusterList(*pNonSeedClusterList);

        for (ClusterList::const_iterator iter = seedClusterList.begin(), iterEnd = seedClusterList.end(); iter != iterEnd; ++iter)
        {
            nonSeedClusterList.erase(*iter);
        }

        if (!nonSeedClusterList.empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListName, nonSeedClusterList));
        }
    }

    if (!seedClusterList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_seedClusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedFindingBaseAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLayerSpan(pCluster) < m_minClusterLayers)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedFindingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    float minClusterLength = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    m_minClusterLayers = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterLayers", m_minClusterLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
