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

#include "TH3F.h"
#include "TFile.h"

using namespace pandora;

namespace lar_content
{

CheatingVertexCreationAlgorithm::CheatingVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_vertexXCorrection(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexCreationAlgorithm::Run()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const VertexList *pVertexList(nullptr); std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    TFile infile(m_pathToSCEFile.c_str(),"READ");

    if(!infile.IsOpen())
    {
        std::cout << "Could not find the space charge effect file !\n";
        return STATUS_CODE_NOT_FOUND;
    }

    TH3F* m_hDx = (TH3F*)infile.Get("hDx");
    TH3F* m_hDy = (TH3F*)infile.Get("hDy");
    TH3F* m_hDz = (TH3F*)infile.Get("hDz");

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
            continue;

        // SCE 2019: First step is to transform position of the neutrino vertex to match SCE maps provided
        const CartesianVector transformedPosition(TransformX(pMCParticle->GetEndpoint().GetX()), TransformY(pMCParticle->GetEndpoint().GetY()),TransformZ(pMCParticle->GetEndpoint().GetZ()) );

        // Then, calculate offsets with the TH3F provided, in the transformed position
        const CartesianVector positionOffset(m_hDx->Interpolate(transformedPosition.GetX(),transformedPosition.GetY(),transformedPosition.GetZ()),
                                             m_hDy->Interpolate(transformedPosition.GetX(),transformedPosition.GetY(),transformedPosition.GetZ()),
                                             m_hDz->Interpolate(transformedPosition.GetX(),transformedPosition.GetY(),transformedPosition.GetZ()));

        // Correct neutrino vertex with the offsets calculated
        const CartesianVector targetVertexSCE(pMCParticle->GetEndpoint().GetX() - positionOffset.GetX(),
                                              pMCParticle->GetEndpoint().GetY() + positionOffset.GetY(),
                                              pMCParticle->GetEndpoint().GetZ() + positionOffset.GetZ());

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = CartesianVector(targetVertexSCE.GetX() + m_vertexXCorrection, targetVertexSCE.GetY(), targetVertexSCE.GetZ());
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

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
/// Transformation coordinates for the maps provided for MicroBooNE as of August 2019

float CheatingVertexCreationAlgorithm::TransformX(float xPosition) const
{
  const float newPos(2.50f - (2.50f/2.56f)*(xPosition/100.0f));
  return newPos;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/// Transform Y to SCE Y coordinate:  [0.0,6.08] --> [0.0,6.0]

float CheatingVertexCreationAlgorithm::TransformY(float yPosition) const
{
    const float newPos((2.50f/2.33f)*((yPosition/100.0f)+1.165f));
    return newPos;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/// Transform Z to SCE Z coordinate:  [0.0,6.97] --> [0.0,7.2]

float CheatingVertexCreationAlgorithm::TransformZ(float zPosition) const
{
    const float newPos((10.0f/10.37f)*(zPosition/100.0f));
    return newPos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "PathToSCEFile", m_pathToSCEFile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexXCorrection", m_vertexXCorrection));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
