/**
 *  @file   LArContent/src/LArTwoDReco/ClusterSplitting/ClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ClusterSplittingAlgorithm::Run()
{
    std::string originalListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalListName));

    if (!m_inputClusterList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_inputClusterList));

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    ClusterSplittingList internalClusterList(pClusterList->begin(), pClusterList->end());

    for (ClusterSplittingList::iterator iter = internalClusterList.begin(); iter != internalClusterList.end(); ++iter)
    {
        Cluster *pCluster = *iter;
        ClusterSplittingList clusterSplittingList;

        if (STATUS_CODE_SUCCESS != this->SplitCluster(pCluster, clusterSplittingList))
            continue;

        internalClusterList.splice(internalClusterList.end(), clusterSplittingList);
        *iter = NULL;
    }

    if (!m_inputClusterList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::SplitCluster(Cluster *const pCluster, ClusterSplittingList &clusterSplittingList) const
{
    // Split cluster into two CaloHit lists
    PandoraContentApi::Cluster::Parameters firstParameters, secondParameters;

    if (STATUS_CODE_SUCCESS != this->DivideCaloHits(pCluster, firstParameters.m_caloHitList, secondParameters.m_caloHitList))
        return STATUS_CODE_NOT_FOUND;

    if (firstParameters.m_caloHitList.empty() || secondParameters.m_caloHitList.empty())
        return STATUS_CODE_NOT_ALLOWED;

    // Begin cluster fragmentation operations
    ClusterList clusterList;
    clusterList.insert(pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    Cluster *pFirstCluster(NULL), *pSecondCluster(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, firstParameters, pFirstCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, secondParameters, pSecondCluster));

    clusterSplittingList.push_back(pFirstCluster);
    clusterSplittingList.push_back(pSecondCluster);

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputClusterList", m_inputClusterList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
