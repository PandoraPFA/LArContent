/**
 *  @file   larpandoracontent/LArUtility/ListDeletionAlgorithm.cc
 *
 *  @brief  Implementation of the list deletion algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ListDeletionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ListDeletionAlgorithm::Run()
{
    for (const std::string &listName : m_pfoListNames)
    {
        const PfoList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pList));

        if (pList && !pList->empty())
        {
            const PfoList listCopy(*pList);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &listCopy, listName));
        }
    }

    for (const std::string &listName : m_clusterListNames)
    {
        const ClusterList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pList));

        if (pList && !pList->empty())
        {
            const ClusterList listCopy(*pList);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &listCopy, listName));
        }
    }

    for (const std::string &listName : m_vertexListNames)
    {
        const VertexList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pList));

        if (pList && !pList->empty())
        {
            const VertexList listCopy(*pList);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &listCopy, listName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListDeletionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "VertexListNames", m_vertexListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
