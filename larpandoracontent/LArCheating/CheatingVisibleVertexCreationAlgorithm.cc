/**
 *  @file   larpandoracontent/LArCheating/CheatingVisibleVertexCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingVisibleVertexCreationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingVisibleVertexCreationAlgorithm::CheatingVisibleVertexCreationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVisibleVertexCreationAlgorithm::Run()
{
    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CheatingVisibleVertexCreationAlgorithm - Could not access pfo list with name " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
            this->SetVertex(pPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingVisibleVertexCreationAlgorithm::SetVertex(const ParticleFlowObject *const pPfo) const
{
    if (!pPfo->GetVertexList().empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const MCParticle *pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));

    const CartesianVector mcVertexPosition = pMCParticle->GetVertex();
    CartesianVector closestPosition(mcVertexPosition);

    CaloHitList spacePoints;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, spacePoints);

    float closestDistance(std::numeric_limits<float>::max());
                    
    for (const CaloHit *const spacePoint : spacePoints)
    {
        float distance((mcVertexPosition - spacePoint->GetPositionVector()).GetMagnitude());

        if (distance < closestDistance)
        {
            closestDistance = distance;
            closestPosition = spacePoint->GetPositionVector();
        }
    }

    const VertexList *pVertexList = NULL;
    std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = closestPosition;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pPfo, pVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVisibleVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
