/**
 *  @file   larpandoracontent/LArUtility/ListMovingAlgorithm.cc
 * 
 *  @brief  Implementation of the list moving algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ListMovingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ListMovingAlgorithm::Run()
{
    for (const std::string &listName : m_pfoListNames)
    {
        const PfoList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pList));

        if (pList && !pList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, listName, m_prefix + listName));
        }
        else if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "ListMovingAlgorithm: unable to find pfos in list " << listName << std::endl;
        }
    }

    for (const std::string &listName : m_clusterListNames)
    {
        const ClusterList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pList));

        if (pList && !pList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, listName, m_prefix + listName));
        }
        else if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "ListMovingAlgorithm: unable to find clusters in list " << listName << std::endl;
        }
    }

    for (const std::string &listName : m_vertexListNames)
    {
        const VertexList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pList));

        if (pList && !pList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, listName, m_prefix + listName));
        }
        else if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "ListMovingAlgorithm: unable to find vertices in list " << listName << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListMovingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "VertexListNames", m_vertexListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "Prefix", m_prefix));

    if (m_prefix.empty())
        return STATUS_CODE_INVALID_PARAMETER;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
