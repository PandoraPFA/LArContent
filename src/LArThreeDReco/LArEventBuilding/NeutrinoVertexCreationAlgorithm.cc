/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/NeutrinoVertexCreationAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/NeutrinoVertexCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoVertexCreationAlgorithm::NeutrinoVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoVertexCreationAlgorithm::Run()
{
    try
    {
        VertexList inputVertexList;
        this->GetInputVertexList(inputVertexList);

        PfoList daughterPfoList;
        this->GetDaughterPfoList(daughterPfoList);

        // TODO: Determine 3D neutrino vertex positions based on daughter Pfos and input vertex
        VertexList outputVertexList(inputVertexList.begin(), inputVertexList.end());

        if (!outputVertexList.empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, outputVertexList));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoVertexCreationAlgorithm::GetInputVertexList(VertexList &vertexList) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputVertexListName, pVertexList));

    if (NULL != pVertexList)
    {
        vertexList.insert(pVertexList->begin(), pVertexList->end());
    }
    else
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoVertexCreationAlgorithm: vertex list " << m_inputVertexListName << " unavailable." << std::endl;
    }

    if (vertexList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoVertexCreationAlgorithm::GetDaughterPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_daughterPfoListNames.begin(), iterEnd = m_daughterPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList(NULL);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, *iter, pPfoList))
        {
            pfoList.insert(pPfoList->begin(), pPfoList->end());
        }
        else
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "NeutrinoVertexCreationAlgorithm: pfo list " << *iter << " unavailable." << std::endl;
        }
    }

    if (pfoList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterPfoListNames", m_daughterPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
