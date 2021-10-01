/**
 *  @file   larpandoracontent/LArCheating/CheatingVertexCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArCheating/CheatingVertexCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingVertexCreationAlgorithm::CheatingVertexCreationAlgorithm() : m_replaceCurrentVertexList(true), m_vertexXCorrection(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexCreationAlgorithm::Run()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    bool cheatVertex(false);
    for (const int nuPdg : m_nuToCheatList)
    {
        if (LArMCParticleHelper::IsCCNuEvent(pMCParticleList, nuPdg))
            cheatVertex = true;
    }

    if (!cheatVertex && !m_nuToCheatList.empty())
        return STATUS_CODE_SUCCESS;

    const VertexList *pVertexList(nullptr);
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
            continue;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = CartesianVector(
            pMCParticle->GetEndpoint().GetX() + m_vertexXCorrection, pMCParticle->GetEndpoint().GetY(), pMCParticle->GetEndpoint().GetZ());
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        CartesianVector position( pMCParticle->GetEndpoint().GetX() + m_vertexXCorrection, pMCParticle->GetEndpoint().GetY(), pMCParticle->GetEndpoint().GetZ());
        CartesianVector origin(0,0,0);

        const Vertex *pVertex(nullptr);
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "NuToCheatList", m_nuToCheatList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexXCorrection", m_vertexXCorrection));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
