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

//------------------------------------------------------------------------------------------------------------------------------------------

TestBeamParticleCreationAlgorithm::TestBeamParticleCreationAlgorithm() :
    m_parentPfoListName(""),
    m_parentVertexListName(""),
    m_keepInteractionVertex(false),
    m_keepStartVertex(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_parentPfoListName, pPfoList));

    PfoList neutrinoPfos;

    for (const Pfo *const pNuPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pNuPfo))
            continue;

        // ATTN: If the test beam score is not set in the neutrino pfo, set the score to 1 as it has been reconstructed under the neutrino hypothesis
        const PropertiesMap properties(pNuPfo->GetPropertiesMap());
        PropertiesMap::const_iterator propertiesIter(properties.find("TestBeamScore"));
        const float testBeamScore(propertiesIter != properties.end() ? (*propertiesIter).second : 1.f);

        const PfoList &nuDaughterPfoList(pNuPfo->GetDaughterPfoList());
        const Pfo *pTestBeamPfo(nullptr);
        CartesianVector positionMinZCaloHit(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

        for (const Pfo *const pNuDaughterPfo : nuDaughterPfoList)
        {
            CaloHitList collectedHits;
            LArPfoHelper::GetCaloHits(pNuDaughterPfo, TPC_3D, collectedHits);

            for (const CaloHit *const pCaloHit : collectedHits)
            {
                if (pCaloHit->GetPositionVector().GetZ() < positionMinZCaloHit.GetZ())
                {
                    positionMinZCaloHit = pCaloHit->GetPositionVector();
                    pTestBeamPfo = pNuDaughterPfo;
                }
            }
        }

        for (const Pfo *const pNuDaughterPfo : nuDaughterPfoList)
        {
            if (pNuDaughterPfo == pTestBeamPfo)
                continue;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pTestBeamPfo, pNuDaughterPfo));
        }

        PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
        pfoMetadata.m_propertiesToAdd["IsTestBeam"] = 1.f;
        pfoMetadata.m_propertiesToAdd["TestBeamScore"] = testBeamScore;

        // ATTN: If the primary pfo is shower like, the target beam particle is most likely an electron/positron.  If the primary pfo is track like, the target
        // beam particle is most likely a pion as pion interactions are more frequent than proton, kaon and muon interactions in the CERN test beam.
        if (std::abs(pTestBeamPfo->GetParticleId()) != E_MINUS)
            pfoMetadata.m_particleId = PI_PLUS;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pTestBeamPfo, pfoMetadata));

        const std::string &daughterPfoListName(LArPfoHelper::IsTrack(pTestBeamPfo) ? m_trackPfoListName : m_showerPfoListName);
        PfoList pfoList(1, pTestBeamPfo);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, daughterPfoListName, m_parentPfoListName, pfoList));

        try
        {
            const Vertex *const pNuDaughterVertex(LArPfoHelper::GetVertex(pTestBeamPfo));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pTestBeamPfo, pNuDaughterVertex));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pNuDaughterVertex, m_daughterVertexListName));

            std::string vertexListName;
            const VertexList *pVertexList(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

            if (m_keepStartVertex)
            {
                PandoraContentApi::Vertex::Parameters parameters;
                parameters.m_position = positionMinZCaloHit;
                parameters.m_vertexLabel = VERTEX_START;
                parameters.m_vertexType = VERTEX_3D;

                const Vertex *pTestBeamStartVertex(nullptr);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pTestBeamStartVertex));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pTestBeamPfo, pTestBeamStartVertex));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_parentVertexListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_parentVertexListName));
            }
            else
            {
                const Vertex *const pNuVertex(LArPfoHelper::GetVertex(pNuPfo));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pNuPfo, pNuVertex));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pTestBeamPfo, pNuVertex));
            }
        }
        catch (const StatusCodeException &)
        {
            std::cout << "TestBeamParticleCreationAlgorithm::Run - unable to modify test beam particle vertex" << std::endl;
        }

        neutrinoPfos.push_back(pNuPfo);
    }

    for (const Pfo *const pNuPfo : neutrinoPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pNuPfo, m_parentPfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentPfoListName", m_parentPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KeepInteractionVertex", m_keepInteractionVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KeepStartVertex", m_keepStartVertex));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentVertexListName", m_parentVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterVertexListName", m_daughterVertexListName));

    if (m_keepInteractionVertex == m_keepStartVertex)
    {
        std::cout << "TestBeamParticleCreationAlgorithm::ReadSettings - must persist one vertex per pfo" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
