/**
 *  @file   larpandoracontent/LArCustomParticles/PCAShowerParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D shower building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPCAHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArShowerPfo.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArCustomParticles/PCAShowerParticleBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PCAShowerParticleBuildingAlgorithm::PCAShowerParticleBuildingAlgorithm() :
    m_cosmicMode(false),
    m_layerFitHalfWindow(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PCAShowerParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const
{
    try
    {
        // In cosmic mode, build showers from all daughter pfos, otherwise require that pfo is shower-like
        if (m_cosmicMode)
        {
            if (LArPfoHelper::IsFinalState(pInputPfo))
                return;
        }
        else
        {
            if (!LArPfoHelper::IsShower(pInputPfo))
                return;
        }

        // Need an input vertex to provide a shower propagation direction
        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        // Need the 3D clusters and hits to calculate PCA components
        ClusterList threeDClusterList;
        LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

        if (threeDClusterList.empty())
            return;

        CaloHitList threeDCaloHitList; // Might be able to get through LArPfoHelper instead
        const Cluster *const pThreeDCluster(threeDClusterList.front());
        pThreeDCluster->GetOrderedCaloHitList().FillCaloHitList(threeDCaloHitList);

        if (threeDCaloHitList.empty())
            return;

        // Run the PCA analysis
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPCAHelper::EigenVectors eigenVecs;
        LArPCAHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPCAHelper::RunPCA(threeDCaloHitList, centroid, eigenValues, eigenVecs);

        // Build a new pfo
        LArShowerPfoFactory pfoFactory;
        LArShowerPfoParameters pfoParameters;
        pfoParameters.m_particleId = (LArPfoHelper::IsShower(pInputPfo) ? pInputPfo->GetParticleId() : E_MINUS);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = pInputPfo->GetMomentum();
        pfoParameters.m_showerVertex = pInputVertex->GetPosition();
        pfoParameters.m_showerCentroid = centroid;

        // Convention: principal axis always points away from vertex (which is assumed to be closer to start of shower than centroid)
        const float testProjection(eigenVecs.at(0).GetDotProduct(pfoParameters.m_showerVertex.Get() - pfoParameters.m_showerCentroid.Get()));
        const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

        pfoParameters.m_showerDirection = eigenVecs.at(0) * directionScaleFactor;
        pfoParameters.m_showerSecondaryVector = eigenVecs.at(1);
        pfoParameters.m_showerTertiaryVector = eigenVecs.at(2);
        pfoParameters.m_showerEigenValues = eigenValues;
        pfoParameters.m_showerLength = this->ShowerLength(pfoParameters.m_showerEigenValues.Get());
        pfoParameters.m_showerOpeningAngle = this->OpeningAngle(pfoParameters.m_showerDirection.Get(),
            pfoParameters.m_showerSecondaryVector.Get(), pfoParameters.m_showerEigenValues.Get());

        const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const ThreeDSlidingFitResult threeDFitResult(pThreeDCluster, m_layerFitHalfWindow, layerPitch);
        pfoParameters.m_showerMinLayerPosition = threeDFitResult.GetGlobalMinLayerPosition();
        pfoParameters.m_showerMaxLayerPosition = threeDFitResult.GetGlobalMaxLayerPosition();

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo,
            pfoFactory));

        const LArShowerPfo *const pLArPfo = dynamic_cast<const LArShowerPfo*>(pOutputPfo);

        if (!pLArPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Build a new vertex
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

CartesianVector PCAShowerParticleBuildingAlgorithm::ShowerLength(const CartesianVector &eigenValues) const
{
    float sl[3] = {0.f, 0.f, 0.f};

    if (eigenValues.GetX() > 0.f)
    {
        sl[0] = 6.f * std::sqrt(eigenValues.GetX());
    }
    else
    {
        std::cout << "The principal eigenvalue is equal to or less than 0." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }

    if (eigenValues.GetY() > 0.f)
        sl[1] = 6.f * std::sqrt(eigenValues.GetY());

    if (eigenValues.GetZ() > 0.f)
        sl[2] = 6.f * std::sqrt(eigenValues.GetZ());

    return CartesianVector(sl[0], sl[1], sl[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PCAShowerParticleBuildingAlgorithm::OpeningAngle(const CartesianVector &principal, const CartesianVector &secondary,
    const CartesianVector &eigenValues) const
{
    const float principalMagnitude(principal.GetMagnitude());
    const float secondaryMagnitude(secondary.GetMagnitude());

    if (std::fabs(principalMagnitude) < std::numeric_limits<float>::epsilon())
    {
        std::cout << "PCAShowerParticleBuildingAlgorithm::OpeningAngle - The principal eigenvector is 0." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
    else if (std::fabs(secondaryMagnitude) < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    const float cosTheta(principal.GetDotProduct(secondary) / (principalMagnitude * secondaryMagnitude));

    if (cosTheta > 1.f)
    {
        std::cout << "PCAShowerParticleBuildingAlgorithm::OpeningAngle - cos(theta) reportedly greater than 1." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const float sinTheta(std::sqrt(1.f - cosTheta * cosTheta));

    if (std::fabs(eigenValues.GetX()) < std::numeric_limits<float>::epsilon())
    {
        std::cout << "PCAShowerParticleBuildingAlgorithm::OpeningAngle - principal eigenvalue less than or equal to 0." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }
    else if (std::fabs(eigenValues.GetY()) < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    return std::atan(std::sqrt(eigenValues.GetY()) * sinTheta / std::sqrt(eigenValues.GetX()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PCAShowerParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosmicMode", m_cosmicMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
