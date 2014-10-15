/**
 *  @file   LArContent/src/LArCheating/CheatingVertexCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating vertex creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArMCParticleHelper.h"

#include "LArCheating/CheatingVertexCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingVertexCreationAlgorithm::CheatingVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexCreationAlgorithm::Run()
{
    const MCParticleList *pMCParticleList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const VertexList *pVertexList(NULL); std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *pMCParticle(*iter);

        if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
            continue;

        if (!pMCParticle->GetParentList().empty())
            continue;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = pMCParticle->GetEndpoint();
        parameters.m_vertexType = VERTEX_3D;

        Vertex *pVertex(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
