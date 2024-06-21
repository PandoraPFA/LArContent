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
    m_warnIfObjectsUnavailable(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPruningAlgorithm::Run()
{
    for (const std::string &listName : m_pfoListNames)
    {
        try
        {
            const PfoList *pPfoList(nullptr);
            const StatusCode statusCode(PandoraContentApi::GetList(*this, listName, pPfoList));

            if (STATUS_CODE_SUCCESS != statusCode)
                throw StatusCodeException(statusCode);

            const PfoList pfoList(*pPfoList);

            for (const ParticleFlowObject *const pPfo : pfoList)
            {
                if (STATUS_CODE_SUCCESS != PandoraContentApi::Delete(*this, pPfo, listName))
                    std::cout << "ListPruningAlgorithm: Could not delete Pfo." << std::endl;
            }
        }
        catch (StatusCodeException &)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPruningAlgorithm: pfo list " << listName << " unavailable." << std::endl;
        }
    }

    for (const std::string &listName : m_clusterListNames)
    {
        try
        {
            const ClusterList *pClusterList(nullptr);
            const StatusCode statusCode(PandoraContentApi::GetList(*this, listName, pClusterList));

            if (STATUS_CODE_SUCCESS != statusCode)
                throw StatusCodeException(statusCode);

            const ClusterList clusterList(*pClusterList);

            for (const Cluster *const pCluster : clusterList)
            {
                if (!m_warnIfObjectsUnavailable && !pCluster->IsAvailable())
                    continue;

                if (STATUS_CODE_SUCCESS != PandoraContentApi::Delete(*this, pCluster, listName) && m_warnIfObjectsUnavailable)
                    std::cout << "ListPruningAlgorithm: Could not delete Cluster." << std::endl;
            }
        }
        catch (StatusCodeException &)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPruningAlgorithm: cluster list " << listName << " unavailable." << std::endl;
        }
    }

    for (const std::string &listName : m_vertexListNames)
    {
        try
        {
            const VertexList *pVertexList(nullptr);
            const StatusCode statusCode(PandoraContentApi::GetList(*this, listName, pVertexList));

            if (STATUS_CODE_SUCCESS != statusCode)
                throw StatusCodeException(statusCode);

            const VertexList vertexList(*pVertexList);

            for (const Vertex *const pVertex : vertexList)
            {
                if (!m_warnIfObjectsUnavailable && !pVertex->IsAvailable())
                    continue;

                if (STATUS_CODE_SUCCESS != PandoraContentApi::Delete(*this, pVertex, listName) && m_warnIfObjectsUnavailable)
                    std::cout << "ListPruningAlgorithm: Could not delete Vertex." << std::endl;
            }
        }
        catch (StatusCodeException &)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPruningAlgorithm: vertex list " << listName << " unavailable." << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPruningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "VertexListNames", m_vertexListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "WarnIfObjectsUnavailable", m_warnIfObjectsUnavailable));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
