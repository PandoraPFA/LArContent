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

        // ATTN: If the primary pfo is shower like, the target beam particle is most likely an electron/positron.  If the primary pfo is track like, the target
        // beam particle is most likely a pion as pion interactions are more frequent than proton, kaon and muon interactions in the CERN test beam.
        if (std::abs(pPrimaryPfo->GetParticleId()) != E_MINUS)
        {
            PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
            pfoMetadata.m_particleId = PI_PLUS;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPrimaryPfo, pfoMetadata));
        }

        try
        {
            if (m_keepStartVertex)
            {
                if (!m_keepInteractionVertex)
                {
                    const Vertex *const pVertex(LArPfoHelper::GetVertex(pPrimaryPfo));
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPrimaryPfo, pVertex));
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pVertex));
                }

                std::string vertexListName;
                const VertexList *pVertexList(nullptr);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

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
        }
        catch (const StatusCodeException &)
        {
            std::cout << "TestBeamParticleCreationAlgorithm::Run - unable to modify test beam particle vertex" << std::endl;
        }

        neutrinoPfos.push_back(pPfo);
    }

    for (const Pfo *const pNeutrinoPfo : neutrinoPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pNeutrinoPfo));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KeepInteractionVertex", m_keepInteractionVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KeepStartVertex", m_keepStartVertex));

    if (m_keepStartVertex)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    if (m_keepInteractionVertex == m_keepStartVertex)
    {
        std::cout << "TestBeamParticleCreationAlgorithm::ReadSettings - must persist one vertex per pfo" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
