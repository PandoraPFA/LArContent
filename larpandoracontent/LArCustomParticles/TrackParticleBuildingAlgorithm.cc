/**
 *  @file   larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D track building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTrackPfo.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackParticleBuildingAlgorithm::TrackParticleBuildingAlgorithm() :
    m_cosmicMode(false),
    m_slidingFitHalfWindow(20)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const
{
    try
    {
        // Need an input vertex to provide a track propagation direction
        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        // In cosmic mode, build tracks from all parent pfos, otherwise require that pfo is track-like
        if (m_cosmicMode)
        {
            if(!LArPfoHelper::IsFinalState(pInputPfo))
                return;
        }
        else
        {
            if (!LArPfoHelper::IsTrack(pInputPfo))
                return;
        }

        // Calculate sliding fit trajectory
        LArTrackStateVector trackStateVector;
        this->GetSlidingFitTrajectory(pInputPfo, pInputVertex, trackStateVector);

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
        pfoParameters.m_trackStateVector = trackStateVector;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo,
            trackFactory));

        const LArTrackPfo *const pLArPfo = dynamic_cast<const LArTrackPfo*>(pOutputPfo);
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

void TrackParticleBuildingAlgorithm::GetSlidingFitTrajectory(const ParticleFlowObject *const pPfo, const Vertex *const pVertex,
    LArTrackStateVector &trackStateVector) const
{
    // Get 3D clusters (normally there should only be one)
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, clusterList);

    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Get seed direction from 3D clusters (bail out for single-hit clusters)
    CartesianVector minPosition(0.f, 0.f, 0.f), maxPosition(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(clusterList, minPosition, maxPosition);

    if ((maxPosition - minPosition).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CartesianVector seedDirection((maxPosition - minPosition).GetUnitVector());
    const CartesianVector seedPosition((maxPosition + minPosition) * 0.5f);

    const bool isForward((seedDirection.GetDotProduct(seedPosition - pVertex->GetPosition()) > 0.f) ? true : false);
    const float scaleFactor(isForward ? +1.f : -1.f);

    // Calculate trajectory from 3D sliding fits
    LArTrackTrajectory trackTrajectory;
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const unsigned int layerWindow(m_slidingFitHalfWindow);

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        try
        {
            const ThreeDSlidingFitResult slidingFitResult(pCluster, layerWindow, layerPitch);

            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            for (CaloHitList::const_iterator hIter = caloHitList.begin(), hIterEnd = caloHitList.end(); hIter != hIterEnd; ++hIter)
            {
                const CaloHit *const pCaloHit3D = *hIter;
                const CaloHit *const pCaloHit2D = static_cast<const CaloHit*>(pCaloHit3D->GetParentAddress());

                try
                {
                    const float rL(slidingFitResult.GetLongitudinalDisplacement(pCaloHit3D->GetPositionVector()));

                    CartesianVector position(0.f, 0.f, 0.f);
                    const StatusCode positionStatusCode(slidingFitResult.GetGlobalFitPosition(rL, position));

                    if (positionStatusCode != STATUS_CODE_SUCCESS)
                        throw StatusCodeException(positionStatusCode);

                    CartesianVector direction(0.f, 0.f, 0.f);
                    const StatusCode directionStatusCode(slidingFitResult.GetGlobalFitDirection(rL, direction));

                    if (directionStatusCode != STATUS_CODE_SUCCESS)
                        throw StatusCodeException(directionStatusCode);

                    const float projection(seedDirection.GetDotProduct(position - seedPosition));

                    trackTrajectory.push_back(LArTrackTrajectoryPoint(projection * scaleFactor,
                        LArTrackState(position, direction * scaleFactor, pCaloHit2D)));
                }
                catch (StatusCodeException &statusCodeException1)
                {
                    if (STATUS_CODE_FAILURE == statusCodeException1.GetStatusCode())
                        throw statusCodeException1;
                }
            }
        }
        catch (StatusCodeException &statusCodeException2)
        {
            if (STATUS_CODE_FAILURE == statusCodeException2.GetStatusCode())
                throw statusCodeException2;
        }
    }

    // Sort trajectory points by distance along track
    std::sort(trackTrajectory.begin(), trackTrajectory.end(), TrackParticleBuildingAlgorithm::SortByHitProjection);

    for (LArTrackTrajectory::const_iterator tIter = trackTrajectory.begin(), tIterEnd = trackTrajectory.end(); tIter != tIterEnd; ++tIter)
    {
        trackStateVector.push_back(tIter->second);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackParticleBuildingAlgorithm::SortByHitProjection(const LArTrackTrajectoryPoint &lhs, const LArTrackTrajectoryPoint &rhs)
{
    if (lhs.first != rhs.first)
        return (lhs.first < rhs.first);

    if (lhs.second.GetCaloHit() && rhs.second.GetCaloHit())
        return (lhs.second.GetCaloHit()->GetInputEnergy() > rhs.second.GetCaloHit()->GetInputEnergy());

    const float dx(lhs.second.GetPosition().GetX() - rhs.second.GetPosition().GetX());
    const float dy(lhs.second.GetPosition().GetY() - rhs.second.GetPosition().GetY());
    const float dz(lhs.second.GetPosition().GetZ() - rhs.second.GetPosition().GetZ());
    return (dx + dy + dz > 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosmicMode", m_cosmicMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
