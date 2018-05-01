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
    m_vertexLowZ(false)
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
        CartesianVector positionMinZCaloHit(0.f, 0.f, 0.f);

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

        if (m_vertexLowZ)
        {
            for (const Vertex *const pVertex : pPrimaryPfo->GetVertexList())
            {
                pPrimaryPfo->RemoveFromPfo(pVertex);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pVertex);
            }

            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = positionMinZCaloHit;
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;

            const Vertex *pVertex(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
            pPrimaryPfo->AddToPfo(pVertex);
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexAtLowZ", m_vertexLowZ));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
