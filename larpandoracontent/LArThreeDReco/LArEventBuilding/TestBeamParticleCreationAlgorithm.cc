/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/TestBeamParticleCreationAlgorithm.cc
 *
 *  @brief  Implementation of the test beam particle creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/TestBeamParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode TestBeamParticleCreationAlgorithm::Run()
{
    const PfoList *pParentNuPfoList(nullptr);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_parentPfoListName, pParentNuPfoList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TestBeamParticleCreationAlgorithm: pfo list " << m_parentPfoListName << " unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    PfoList neutrinoPfos;

    for (const Pfo *const pNuPfo : *pParentNuPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pNuPfo))
            continue;

        neutrinoPfos.push_back(pNuPfo);

        const Pfo *pTestBeamPfo(nullptr);
        CartesianVector testBeamStartVertex(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SetupTestBeamPfo(pNuPfo, pTestBeamPfo, testBeamStartVertex));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SetupTestBeamVertex(pNuPfo, pTestBeamPfo, testBeamStartVertex));
    }

    for (const Pfo *const pNuPfo : neutrinoPfos)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pNuPfo, m_parentPfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::SetupTestBeamPfo(const Pfo *const pNuPfo, const Pfo *&pTestBeamPfo, CartesianVector &testBeamStartVertex) const
{
    pTestBeamPfo = nullptr;
    testBeamStartVertex = CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    for (const Pfo *const pNuDaughterPfo : pNuPfo->GetDaughterPfoList())
    {
        CaloHitList collectedHits;
        LArPfoHelper::GetCaloHits(pNuDaughterPfo, TPC_3D, collectedHits);

        for (const CaloHit *const pCaloHit : collectedHits)
        {
            if (pCaloHit->GetPositionVector().GetZ() < testBeamStartVertex.GetZ())
            {
                testBeamStartVertex = pCaloHit->GetPositionVector();
                pTestBeamPfo = pNuDaughterPfo;
            }
        }
    }

    if (!pTestBeamPfo)
        return STATUS_CODE_NOT_FOUND;

    for (const Pfo *const pNuDaughterPfo : pNuPfo->GetDaughterPfoList())
    {
        if (pNuDaughterPfo != pTestBeamPfo)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pTestBeamPfo, pNuDaughterPfo));
    }

    // Move test beam pfo to parent list from its initial track or shower list
    const std::string &originalListName(LArPfoHelper::IsTrack(pTestBeamPfo) ? m_trackPfoListName : m_showerPfoListName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, originalListName, m_parentPfoListName, PfoList(1, pTestBeamPfo)));

    // Alter metadata: test beam score (1, if not previously set) and particle id (electron if primary pfo shower-like, else charged pion)
    PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
    pfoMetadata.m_propertiesToAdd["IsTestBeam"] = 1.f;
    pfoMetadata.m_propertiesToAdd["TestBeamScore"] = (pNuPfo->GetPropertiesMap().count("TestBeamScore")) ? pNuPfo->GetPropertiesMap().at("TestBeamScore") : 1.f;
    pfoMetadata.m_particleId = (std::fabs(pTestBeamPfo->GetParticleId()) == E_MINUS) ? E_MINUS : PI_PLUS;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pTestBeamPfo, pfoMetadata));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::SetupTestBeamVertex(const Pfo *const pNuPfo, const Pfo *const pTestBeamPfo, const CartesianVector &testBeamStartVertex) const
{
    try
    {
        const Vertex *const pInitialVertex(LArPfoHelper::GetVertex(pTestBeamPfo));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pTestBeamPfo, pInitialVertex));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pInitialVertex, m_daughterVertexListName));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "TestBeamParticleCreationAlgorithm::SetupTestBeamVertex - Test beam particle has no initial vertex" << std::endl;
    }

    // Set start vertex
    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = testBeamStartVertex;
    parameters.m_vertexLabel = VERTEX_START;
    parameters.m_vertexType = VERTEX_3D;
    const Vertex *pTestBeamStartVertex(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pTestBeamStartVertex));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_parentVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pTestBeamPfo, pTestBeamStartVertex));

    // Retain interaction vertex
    try
    {
        const Vertex *const pNuVertex(LArPfoHelper::GetVertex(pNuPfo));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pNuPfo, pNuVertex));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pTestBeamPfo, pNuVertex));
        // ATTN This vertex already lives in m_parentVertexListName
    }
    catch (const StatusCodeException &)
    {
        std::cout << "TestBeamParticleCreationAlgorithm::SetupTestBeamVertex - Cannot transfer interaction vertex to test beam particle" << std::endl;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_parentVertexListName));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ParentPfoListName", m_parentPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackPfoListName", m_trackPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ParentVertexListName", m_parentVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "DaughterVertexListName", m_daughterVertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
