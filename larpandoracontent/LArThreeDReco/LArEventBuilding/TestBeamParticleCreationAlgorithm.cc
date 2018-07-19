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
    m_pfoListName(""),
    m_vertexListName(""),
    m_keepInteractionVertex(false),
    m_keepStartVertex(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    PfoList neutrinoPfos;

    for (const Pfo *const pPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pPfo))
            continue;

        // ATTN: If the test beam score is not set in the neutrino pfo, set the score to 1 as it has been reconstructed under the neutrino hypothesis
        const PropertiesMap properties(pPfo->GetPropertiesMap());
        PropertiesMap::const_iterator propertiesIter(properties.find("TestBeamScore"));
        const float testBeamScore(propertiesIter != properties.end() ? (*propertiesIter).second : 1.f);

        const PfoList &daughterList(pPfo->GetDaughterPfoList());
        const Pfo *pPrimaryPfo(nullptr);
        CartesianVector positionMinZCaloHit(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

        for (const Pfo *const pDaughterPfo : daughterList)
        {
            CaloHitList collectedHits;
            LArPfoHelper::GetCaloHits(pDaughterPfo, TPC_3D, collectedHits);

            for (const CaloHit *const pCaloHit : collectedHits)
            {
                if (pCaloHit->GetPositionVector().GetZ() < positionMinZCaloHit.GetZ())
                {
                    positionMinZCaloHit = pCaloHit->GetPositionVector();
                    pPrimaryPfo = pDaughterPfo;
                }
            }
        }

        for (const Pfo *const pPrimaryDaughterPfo : daughterList)
        {
            if (pPrimaryDaughterPfo == pPrimaryPfo)
                continue;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPrimaryPfo, pPrimaryDaughterPfo));
        }

        PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
        pfoMetadata.m_propertiesToAdd["IsTestBeam"] = 1.f;
        pfoMetadata.m_propertiesToAdd["TestBeamScore"] = testBeamScore;

        // ATTN: If the primary pfo is shower like, the target beam particle is most likely an electron/positron.  If the primary pfo is track like, the target
        // beam particle is most likely a pion as pion interactions are more frequent than proton, kaon and muon interactions in the CERN test beam.
        if (std::abs(pPrimaryPfo->GetParticleId()) != E_MINUS)
            pfoMetadata.m_particleId = PI_PLUS;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPrimaryPfo, pfoMetadata));

        const std::string previousListName(LArPfoHelper::IsTrack(pPrimaryPfo) ? m_trackPfoListName : m_showerPfoListName);
        PfoList pfoList;
        pfoList.push_back(pPrimaryPfo);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, previousListName, m_pfoListName, pfoList));

        try
        {
            const Vertex *const pOriginalVertex(LArPfoHelper::GetVertex(pPrimaryPfo));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPrimaryPfo, pOriginalVertex));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pOriginalVertex, m_daughterVertexListName));

            std::string vertexListName;
            const VertexList *pVertexList(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

            if (m_keepStartVertex)
            {
                PandoraContentApi::Vertex::Parameters parameters;
                parameters.m_position = positionMinZCaloHit;
                parameters.m_vertexLabel = VERTEX_START;
                parameters.m_vertexType = VERTEX_3D;

                const Vertex *pVertex(nullptr);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPrimaryPfo, pVertex));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_vertexListName));
            }
            else
            {
                const Vertex *const pVertex(LArPfoHelper::GetVertex(pPfo));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pVertex));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPrimaryPfo, pVertex));
            }
        }
        catch (const StatusCodeException &)
        {
            std::cout << "TestBeamParticleCreationAlgorithm::Run - unable to modify test beam particle vertex" << std::endl;
        }

        neutrinoPfos.push_back(pPfo);
    }

    for (const Pfo *const pNeutrinoPfo : neutrinoPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pNeutrinoPfo, m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KeepInteractionVertex", m_keepInteractionVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KeepStartVertex", m_keepStartVertex));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterVertexListName", m_daughterVertexListName));

    if (m_keepInteractionVertex == m_keepStartVertex)
    {
        std::cout << "TestBeamParticleCreationAlgorithm::ReadSettings - must persist one vertex per pfo" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
