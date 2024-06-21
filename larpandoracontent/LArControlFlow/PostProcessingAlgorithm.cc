/**
 *  @file   larpandoracontent/LArControlFlow/PostProcessingAlgorithm.cc
 *
 *  @brief  Implementation of the list moving algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/PostProcessingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PostProcessingAlgorithm::PostProcessingAlgorithm() :
    m_listCounter(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PostProcessingAlgorithm::Reset()
{
    m_listCounter = 0;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PostProcessingAlgorithm::Run()
{
    for (const std::string &listName : m_pfoListNames)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RenameList<PfoList>(listName));

    for (const std::string &listName : m_clusterListNames)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RenameList<ClusterList>(listName));

    for (const std::string &listName : m_vertexListNames)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RenameList<VertexList>(listName));

    for (const std::string &listName : m_caloHitListNames)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RenameList<CaloHitList>(listName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<CaloHit>(*this));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Track>(*this));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<MCParticle>(*this));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Cluster>(*this));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<ParticleFlowObject>(*this));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Vertex>(*this));

    if (!m_currentPfoListReplacement.empty())
    {
        const std::string replacementListName(m_currentPfoListReplacement + TypeToString(m_listCounter));

        if (STATUS_CODE_SUCCESS != PandoraContentApi::ReplaceCurrentList<Pfo>(*this, replacementListName))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PostProcessingAlgorithm: could not replace current pfo list with list named: " << replacementListName << std::endl;

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Pfo>(*this));
        }
    }

    ++m_listCounter;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode PostProcessingAlgorithm::RenameList(const std::string &oldListName) const
{
    const std::string newListName(oldListName + TypeToString(m_listCounter));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, PandoraContentApi::RenameList<T>(*this, oldListName, newListName));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PostProcessingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "VertexListNames", m_vertexListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CurrentPfoListReplacement", m_currentPfoListReplacement));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template pandora::StatusCode PostProcessingAlgorithm::RenameList<PfoList>(const std::string &) const;
template pandora::StatusCode PostProcessingAlgorithm::RenameList<ClusterList>(const std::string &) const;
template pandora::StatusCode PostProcessingAlgorithm::RenameList<VertexList>(const std::string &) const;
template pandora::StatusCode PostProcessingAlgorithm::RenameList<CaloHitList>(const std::string &) const;

} // namespace lar_content
