/**
 *  @file   LArContent/src/LArUtility/ListChangingAlgorithm.cc
 * 
 *  @brief  Implementation of the list changing algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ListChangingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ListChangingAlgorithm::Run()
{
    if (!m_caloHitListName.empty())
    {
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_caloHitListName));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            return statusCode;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo() && (STATUS_CODE_NOT_FOUND == statusCode))
            std::cout << "ListChangingAlgorithm: calohit list not found " << m_caloHitListName << std::endl;
    }

    if (!m_clusterListName.empty())
    {
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            return statusCode;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo() && (STATUS_CODE_NOT_FOUND == statusCode))
            std::cout << "ListChangingAlgorithm: cluster list not found " << m_clusterListName << std::endl;
    }

    if (!m_vertexListName.empty())
    {
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_vertexListName));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            return statusCode;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo() && (STATUS_CODE_NOT_FOUND == statusCode))
            std::cout << "ListChangingAlgorithm: vertex list not found " << m_vertexListName << std::endl;
    }

    if (!m_pfoListName.empty())
    {
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_pfoListName));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            return statusCode;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo() && (STATUS_CODE_NOT_FOUND == statusCode))
            std::cout << "ListChangingAlgorithm: pfo list not found " << m_pfoListName << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListChangingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
