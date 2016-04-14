/**
 *  @file   LArContent/src/LArTwoDReco/ClusterSplitting/ClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ClusterSplittingAlgorithm::Run()
{
    std::string originalListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalListName));

    if (!m_inputClusterList.empty())
    {
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_inputClusterList));

        if (STATUS_CODE_NOT_FOUND == statusCode)
        {
            std::cout << "ClusterSplittingAlgorithm: cluster list not found " << m_inputClusterList << std::endl;
            return STATUS_CODE_SUCCESS;
        }

        if (STATUS_CODE_SUCCESS != statusCode)
            return statusCode;
    }

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    ClusterSplittingList internalClusterList(pClusterList->begin(), pClusterList->end());
    internalClusterList.sort(LArClusterHelper::SortByNHits);

    for (ClusterSplittingList::iterator iter = internalClusterList.begin(); iter != internalClusterList.end(); ++iter)
    {
        const Cluster *const pCluster = *iter;
        ClusterSplittingList clusterSplittingList;

        if (STATUS_CODE_SUCCESS != this->SplitCluster(pCluster, clusterSplittingList))
            continue;

        internalClusterList.splice(internalClusterList.end(), clusterSplittingList);
        *iter = NULL;
    }

    if (!m_inputClusterList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalListName));

const ClusterList *pClusterList1(nullptr);
if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ClustersU", pClusterList1))
{
    ClusterVector clusterVector1(pClusterList1->begin(), pClusterList1->end());
    std::sort(clusterVector1.begin(), clusterVector1.end(), LArClusterHelper::SortByNHits);
    for (const Cluster *const pCluster1 : clusterVector1)
        std::cout << "Alg " << this->GetType() << "Cluster " << pCluster1->GetNCaloHits() << ", E " << pCluster1->GetHadronicEnergy() << std::endl;
}
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, ClusterSplittingList &clusterSplittingList) const
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
    const Cluster *pFirstCluster(NULL), *pSecondCluster(NULL);
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
