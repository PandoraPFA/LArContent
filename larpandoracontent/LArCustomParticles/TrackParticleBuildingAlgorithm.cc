/**
 *  @file   larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D track building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Managers/GeometryManager.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTrackPfo.h"

#include "larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackParticleBuildingAlgorithm::TrackParticleBuildingAlgorithm() :
    m_slidingFitHalfWindow(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *&pOutputPfo) const
{
    try
    {
        // Need an input vertex to provide a track propagation direction
        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        // In cosmic mode, build tracks from all parent pfos, otherwise require that pfo is track-like
        if (LArPfoHelper::IsNeutrinoFinalState(pInputPfo))
        {
            if (!LArPfoHelper::IsTrack(pInputPfo))
                return;
        }
        else
        {
            if (!LArPfoHelper::IsFinalState(pInputPfo))
                return;

            if (LArPfoHelper::IsNeutrino(pInputPfo))
                return;
        }

        // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
        const LArTPC *const pFirstLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
        const float layerPitch(pFirstLArTPC->GetWirePitchW());

        // Calculate sliding fit trajectory
        LArTrackStateVector trackStateVector;
        LArPfoHelper::GetSlidingFitTrajectory(pInputPfo, pInputVertex, m_slidingFitHalfWindow, layerPitch, trackStateVector);

        if (trackStateVector.empty())
            return;

        // Build track-like pfo from track trajectory (TODO Correct these placeholder parameters)
        LArTrackPfoFactory trackFactory;
        LArTrackPfoParameters pfoParameters;
        pfoParameters.m_particleId = (LArPfoHelper::IsTrack(pInputPfo) ? pInputPfo->GetParticleId() : MU_MINUS);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = pInputPfo->GetMomentum();
        pfoParameters.m_propertiesToAdd = pInputPfo->GetPropertiesMap();
        pfoParameters.m_trackStateVector = std::move(trackStateVector);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo, trackFactory));

        const LArTrackPfo *const pLArPfo = dynamic_cast<const LArTrackPfo *>(pOutputPfo);
        if (NULL == pLArPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Now update vertex and direction
        PandoraContentApi::ParticleFlowObject::Metadata pfodata;
        pfodata.m_momentum = pLArPfo->GetVertexDirection();
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pOutputPfo, pfodata));

        const Vertex *pOutputVertex(NULL);

        PandoraContentApi::Vertex::Parameters vtxParameters;
        vtxParameters.m_position = pLArPfo->GetVertexPosition();
        vtxParameters.m_vertexLabel = pInputVertex->GetVertexLabel();
        vtxParameters.m_vertexType = pInputVertex->GetVertexType();

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, vtxParameters, pOutputVertex));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pOutputPfo, pOutputVertex));
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
