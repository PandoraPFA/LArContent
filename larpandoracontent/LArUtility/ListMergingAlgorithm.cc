/**
 *  @file   larpandoracontent/LArUtility/ListMergingAlgorithm.cc
 *
 *  @brief  Implementation of the list merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ListMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ListMergingAlgorithm::Run()
{
    // Vertex list merging
    if (m_sourceVertexListNames.size() != m_targetVertexListNames.size())
        return STATUS_CODE_FAILURE;

    for (unsigned int iIndex = 0, iIndexEnd = m_sourceVertexListNames.size(); iIndex < iIndexEnd; ++iIndex)
    {
        const std::string &sourceListName(m_sourceVertexListNames.at(iIndex));
        const std::string &targetListName(m_targetVertexListNames.at(iIndex));

        const StatusCode statusCode(PandoraContentApi::SaveList<Vertex>(*this, sourceListName, targetListName));

        if (STATUS_CODE_SUCCESS != statusCode)
        {
            if (STATUS_CODE_NOT_FOUND == statusCode)
            {
                if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListMergingAlgorithm: vertex list not found, source: " << sourceListName
                              << ", target: " << targetListName << std::endl;
            }
            else if (STATUS_CODE_NOT_INITIALIZED == statusCode)
            {
                if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListMergingAlgorithm: no vertexs to move, source: " << sourceListName << ", target: " << targetListName << std::endl;
            }
            else
            {
                std::cout << "ListMergingAlgorithm: error in vertex merging, source: " << sourceListName << ", target: " << targetListName
                          << std::endl;
                return statusCode;
            }
        }
    }

    // Cluster list merging
    if (m_sourceClusterListNames.size() != m_targetClusterListNames.size())
        return STATUS_CODE_FAILURE;

    for (unsigned int iIndex = 0, iIndexEnd = m_sourceClusterListNames.size(); iIndex < iIndexEnd; ++iIndex)
    {
        const std::string &sourceListName(m_sourceClusterListNames.at(iIndex));
        const std::string &targetListName(m_targetClusterListNames.at(iIndex));

        const StatusCode statusCode(PandoraContentApi::SaveList<Cluster>(*this, sourceListName, targetListName));

        if (STATUS_CODE_SUCCESS != statusCode)
        {
            if (STATUS_CODE_NOT_FOUND == statusCode)
            {
                if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListMergingAlgorithm: cluster list not found, source: " << sourceListName
                              << ", target: " << targetListName << std::endl;
            }
            else if (STATUS_CODE_NOT_INITIALIZED == statusCode)
            {
                if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListMergingAlgorithm: no clusters to move, source: " << sourceListName << ", target: " << targetListName << std::endl;
            }
            else
            {
                std::cout << "ListMergingAlgorithm: error in cluster merging, source: " << sourceListName << ", target: " << targetListName
                          << std::endl;
                return statusCode;
            }
        }
    }

    // Pfo list merging
    if (m_sourcePfoListNames.size() != m_targetPfoListNames.size())
        return STATUS_CODE_FAILURE;

    for (unsigned int iIndex = 0, iIndexEnd = m_sourcePfoListNames.size(); iIndex < iIndexEnd; ++iIndex)
    {
        const std::string &sourceListName(m_sourcePfoListNames.at(iIndex));
        const std::string &targetListName(m_targetPfoListNames.at(iIndex));

        const StatusCode statusCode(PandoraContentApi::SaveList<Pfo>(*this, sourceListName, targetListName));

        if (STATUS_CODE_SUCCESS != statusCode)
        {
            if (STATUS_CODE_NOT_FOUND == statusCode)
            {
                if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListMergingAlgorithm: pfo list not found, source: " << sourceListName << ", target: " << targetListName << std::endl;
            }
            else if (STATUS_CODE_NOT_INITIALIZED == statusCode)
            {
                if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListMergingAlgorithm: no pfos to move, source: " << sourceListName << ", target: " << targetListName << std::endl;
            }
            else
            {
                std::cout << "ListMergingAlgorithm: error in pfo merging, source: " << sourceListName << ", target: " << targetListName << std::endl;
                return statusCode;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "SourceVertexListNames", m_sourceVertexListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "TargetVertexListNames", m_targetVertexListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "SourceClusterListNames", m_sourceClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "TargetClusterListNames", m_targetClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "SourcePfoListNames", m_sourcePfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "TargetPfoListNames", m_targetPfoListNames));

    if ((m_sourceClusterListNames.size() != m_targetClusterListNames.size()) || (m_sourcePfoListNames.size() != m_targetPfoListNames.size()))
    {
        std::cout << "ListMergingAlgorithm::ReadSettings: invalid list configuration " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
