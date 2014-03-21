/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/LArSeedFinding/SeedFindingBaseAlgorithm.cc
 *
 *  @brief  Implementation of the seed finding base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArSeedFinding/SeedFindingBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedFindingBaseAlgorithm::Run()
{
    // Organize input lists
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_seedClusterListName, pSeedClusterList));

    const ClusterList *pNonSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    const bool nonSeedListExists(NULL != pNonSeedClusterList);

    if (!nonSeedListExists)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pNonSeedClusterList));

        if (pSeedClusterList == pNonSeedClusterList)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    // Select the new particle seeds
    ClusterVector clusterVector(pNonSeedClusterList->begin(), pNonSeedClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    ClusterList seedClusterList;
    this->GetSeedClusterList(clusterVector, seedClusterList);

    // Organize output lists
    if (!seedClusterList.empty())
    {
        if (nonSeedListExists)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_nonSeedClusterListName, m_seedClusterListName, seedClusterList));
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_seedClusterListName, seedClusterList));
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
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_nonSeedClusterListName, nonSeedClusterList));
        }
    }

    if (!seedClusterList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_seedClusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedFindingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
