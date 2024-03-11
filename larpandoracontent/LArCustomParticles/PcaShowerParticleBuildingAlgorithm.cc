/**
 *  @file   larpandoracontent/LArCustomParticles/PcaShowerParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D shower building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArShowerPfo.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArCustomParticles/PcaShowerParticleBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PcaShowerParticleBuildingAlgorithm::PcaShowerParticleBuildingAlgorithm() :
    m_layerFitHalfWindow(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PcaShowerParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *&pOutputPfo) const
{
    try
    {
        // In cosmic mode, build showers from all daughter pfos, otherwise require that pfo is shower-like
        if (LArPfoHelper::IsNeutrinoFinalState(pInputPfo))
        {
            if (!LArPfoHelper::IsShower(pInputPfo))
                return;
        }
        else
        {
            if (LArPfoHelper::IsFinalState(pInputPfo))
                return;

            if (LArPfoHelper::IsNeutrino(pInputPfo))
                return;
        }

        // Need an input vertex to provide a shower propagation direction
        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        // Run the PCA analysis
        const LArShowerPCA showerPCA(LArPfoHelper::GetPrincipalComponents(pInputPfo, pInputVertex));

        // Build a new pfo
        LArShowerPfoFactory pfoFactory;
        LArShowerPfoParameters pfoParameters;
        pfoParameters.m_particleId = (LArPfoHelper::IsShower(pInputPfo) ? pInputPfo->GetParticleId() : E_MINUS);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = pInputPfo->GetMomentum();
        pfoParameters.m_propertiesToAdd = pInputPfo->GetPropertiesMap();
        pfoParameters.m_showerVertex = pInputVertex->GetPosition();
        pfoParameters.m_showerCentroid = showerPCA.GetCentroid();
        pfoParameters.m_showerDirection = showerPCA.GetPrimaryAxis();
        pfoParameters.m_showerSecondaryVector = showerPCA.GetSecondaryAxis();
        pfoParameters.m_showerTertiaryVector = showerPCA.GetTertiaryAxis();
        pfoParameters.m_showerEigenValues = showerPCA.GetEigenValues();
        pfoParameters.m_showerLength = showerPCA.GetAxisLengths();
        pfoParameters.m_showerOpeningAngle =
            (showerPCA.GetPrimaryLength() > 0.f ? std::atan(showerPCA.GetSecondaryLength() / showerPCA.GetPrimaryLength()) : 0.f);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo, pfoFactory));

        const LArShowerPfo *const pLArPfo = dynamic_cast<const LArShowerPfo *>(pOutputPfo);

        if (!pLArPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Build a new vertex - TODO: tune vertex position based on PCA results
        const Vertex *pOutputVertex = nullptr;

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

StatusCode PcaShowerParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LayerFitHalfWindow", m_layerFitHalfWindow));

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
