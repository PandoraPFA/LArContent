/**
 *  @file   larpandoracontent/LArUtility/ListPruningAlgorithm.cc
 *
 *  @brief  Implementation of the list pruning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ListPruningAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ListPruningAlgorithm::ListPruningAlgorithm() :
    m_warnIfClustersUnavailable(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPruningAlgorithm::Run()
{
    for (StringVector::const_iterator iter = m_pfoListNames.begin(), iterEnd = m_pfoListNames.end(); iter != iterEnd; ++iter)
    {
        const std::string &listName(*iter);

        try
        {
            const PfoList *pPfoList = NULL;
            const StatusCode statusCode(PandoraContentApi::GetList(*this, listName, pPfoList));

            if (STATUS_CODE_SUCCESS != statusCode)
                throw StatusCodeException(statusCode);

            const PfoList pfoList(*pPfoList);

            for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
            {
                if (STATUS_CODE_SUCCESS != PandoraContentApi::Delete(*this, *pIter, listName))
                {
                    std::cout << "ListPruningAlgorithm: Could not delete Pfo." << std::endl;
                }
            }
        }
        catch (StatusCodeException &)
        {
            std::cout << "ListPruningAlgorithm: pfo list " << listName << " unavailable." << std::endl;
        }
    }

    for (StringVector::const_iterator iter = m_clusterListNames.begin(), iterEnd = m_clusterListNames.end(); iter != iterEnd; ++iter)
    {
        const std::string &listName(*iter);

        try
        {
            const ClusterList *pClusterList = NULL;
            const StatusCode statusCode(PandoraContentApi::GetList(*this, listName, pClusterList));

            if (STATUS_CODE_SUCCESS != statusCode)
                throw StatusCodeException(statusCode);

            const ClusterList clusterList(*pClusterList);

            for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
            {
                if (!m_warnIfClustersUnavailable && !(*cIter)->IsAvailable())
                    continue;

                if (STATUS_CODE_SUCCESS != PandoraContentApi::Delete(*this, *cIter, listName))
                {
                    if (m_warnIfClustersUnavailable)
                        std::cout << "ListPruningAlgorithm: Could not delete Pfo." << std::endl;
                }
            }
        }
        catch (StatusCodeException &)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPruningAlgorithm: cluster list " << listName << " unavailable." << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPruningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WarnIfClustersUnavailable", m_warnIfClustersUnavailable));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
