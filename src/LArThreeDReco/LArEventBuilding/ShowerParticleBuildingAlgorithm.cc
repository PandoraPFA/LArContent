/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/ShowerParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D shower building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"
#include "LArObjects/LArShowerPfo.h"

#include "LArThreeDReco/LArEventBuilding/ShowerParticleBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void ShowerParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const
{
    try
    {
        // Need an input vertex to provide a shower propagation direction
        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        // Calculate shower variables if pfo is shower-like
        if (!LArPfoHelper::IsShower(pInputPfo))
            return;

        // Build a new pfo
        LArShowerPfoFactory pfoFactory;
        LArShowerPfoParameters pfoParameters;
        pfoParameters.m_particleId = pInputPfo->GetParticleId();
        pfoParameters.m_charge = pInputPfo->GetCharge();
        pfoParameters.m_mass = pInputPfo->GetMass();
        pfoParameters.m_energy = pInputPfo->GetEnergy();
        pfoParameters.m_momentum = pInputPfo->GetMomentum();
        pfoParameters.m_additionalProperty = "HelloWorld";

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo,
            pfoFactory));

        const LArShowerPfo *const pLArPfo = dynamic_cast<const LArShowerPfo*>(pOutputPfo);
        if (NULL == pLArPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Build a new vertex
        const Vertex *pOutputVertex(NULL);

        PandoraContentApi::Vertex::Parameters vtxParameters;
        vtxParameters.m_position = pInputVertex->GetPosition();
        vtxParameters.m_vertexLabel = pInputVertex->GetVertexLabel();
        vtxParameters.m_vertexType = pInputVertex->GetVertexType();

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, vtxParameters, pOutputVertex));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pOutputPfo, pOutputVertex));
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
