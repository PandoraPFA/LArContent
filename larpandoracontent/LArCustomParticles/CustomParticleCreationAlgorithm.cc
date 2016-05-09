/**
 *  @file   larpandoracontent/LArCustomParticles/CustomParticleCreationAlgorithm.cc
 *
 *  @brief  Implementation of the 3D particle creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCustomParticles/CustomParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CustomParticleCreationAlgorithm::Run()
{
    // Get input Pfo List
    const PfoList *pPfoList(NULL);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_pfoListName, pPfoList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CustomParticleCreationAlgorithm: cannot find pfo list " << m_pfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Get input Vertex List
    const VertexList *pVertexList(NULL);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_vertexListName, pVertexList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CustomParticleCreationAlgorithm: cannot find vertex list " << m_vertexListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Create temporary lists
    const PfoList *pTempPfoList = NULL; std::string tempPfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempPfoList,
        tempPfoListName));

    const VertexList *pTempVertexList = NULL; std::string tempVertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempVertexList,
        tempVertexListName));

    // Loop over input Pfos
    PfoList pfoList(pPfoList->begin(), pPfoList->end());
    VertexList vertexList(pVertexList->begin(), pVertexList->end());

    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pInputPfo = *iter;

        if (pInputPfo->GetVertexList().empty())
            continue;

        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        if (0 == vertexList.count(pInputVertex))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Build a new pfo and vertex from the old pfo
        const ParticleFlowObject *pOutputPfo(NULL);

        this->CreatePfo(pInputPfo, pOutputPfo);

        if (NULL == pOutputPfo)
            continue;

        if (pOutputPfo->GetVertexList().empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Transfer clusters and hierarchy information to new pfo, and delete old pfo and vertex
        ClusterList clusterList(pInputPfo->GetClusterList().begin(), pInputPfo->GetClusterList().end());
        PfoList parentList(pInputPfo->GetParentPfoList().begin(), pInputPfo->GetParentPfoList().end());
        PfoList daughterList(pInputPfo->GetDaughterPfoList().begin(), pInputPfo->GetDaughterPfoList().end());

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pInputPfo, m_pfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pInputVertex, m_vertexListName));

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Cluster>(*this, pOutputPfo, *cIter));
        }

        for (PfoList::const_iterator pIter = parentList.begin(), pIterEnd = parentList.end(); pIter != pIterEnd; ++pIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, *pIter, pOutputPfo));
        }

        for (PfoList::const_iterator dIter = daughterList.begin(), dIterEnd = daughterList.end(); dIter != dIterEnd; ++dIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pOutputPfo, *dIter));
        }
    }

    if (!pTempPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_pfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_pfoListName));
    }

    if (!pTempVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_vertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CustomParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
