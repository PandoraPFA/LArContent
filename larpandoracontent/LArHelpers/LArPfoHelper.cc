/**
 *  @file   larpandoracontent/LArHelpers/LArPfoHelper.cc
 *
 *  @brief  Implementation of the pfo helper class.
 *
 *  $Log: $
 */

#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

void LArPfoHelper::GetCaloHits(const PfoList &pfoList, const HitType &hitType, CaloHitList &caloHitList)
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        LArPfoHelper::GetCaloHits(*pIter, hitType, caloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetCaloHits(const ParticleFlowObject *const pPfo, const HitType &hitType, CaloHitList &caloHitList)
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        (*cIter)->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetIsolatedCaloHits(const ParticleFlowObject *const pPfo, const HitType &hitType, CaloHitList &caloHitList)
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        caloHitList.insert(caloHitList.end(), (*cIter)->GetIsolatedCaloHitList().begin(), (*cIter)->GetIsolatedCaloHitList().end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetClusters(const PfoList &pfoList, const HitType &hitType, ClusterList &clusterList)
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        LArPfoHelper::GetClusters(*pIter, hitType, clusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetClusters(const ParticleFlowObject *const pPfo, const HitType &hitType, ClusterList &clusterList)
{
    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pPfoCluster = *cIter;

        if (hitType != LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        clusterList.push_back(pPfoCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetTwoDClusterList(const ParticleFlowObject *const pPfo, ClusterList &clusterList)
{
    for (ClusterList::const_iterator cIter = pPfo->GetClusterList().begin(), cIterEnd = pPfo->GetClusterList().end();
        cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        clusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetThreeDClusterList(const ParticleFlowObject *const pPfo, ClusterList &clusterList)
{
    for (ClusterList::const_iterator cIter = pPfo->GetClusterList().begin(), cIterEnd = pPfo->GetClusterList().end();
        cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        clusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllConnectedPfos(const PfoList &inputPfoList, PfoList &outputPfoList)
{
    for (PfoList::const_iterator pIter = inputPfoList.begin(), pIterEnd = inputPfoList.end(); pIter != pIterEnd; ++pIter)
    {
        LArPfoHelper::GetAllConnectedPfos(*pIter, outputPfoList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllConnectedPfos(const ParticleFlowObject *const pPfo, PfoList &outputPfoList)
{
    if (outputPfoList.end() != std::find(outputPfoList.begin(), outputPfoList.end(), pPfo))
        return;

    outputPfoList.push_back(pPfo);
    LArPfoHelper::GetAllConnectedPfos(pPfo->GetParentPfoList(), outputPfoList);
    LArPfoHelper::GetAllConnectedPfos(pPfo->GetDaughterPfoList(), outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllDownstreamPfos(const PfoList &inputPfoList, PfoList &outputPfoList)
{
    for (PfoList::const_iterator pIter = inputPfoList.begin(), pIterEnd = inputPfoList.end(); pIter != pIterEnd; ++pIter)
    {
        LArPfoHelper::GetAllDownstreamPfos(*pIter, outputPfoList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetAllDownstreamPfos(const ParticleFlowObject *const pPfo, PfoList &outputPfoList)
{
    if (outputPfoList.end() != std::find(outputPfoList.begin(), outputPfoList.end(), pPfo))
        return;

    outputPfoList.push_back(pPfo);
    LArPfoHelper::GetAllDownstreamPfos(pPfo->GetDaughterPfoList(), outputPfoList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetTwoDLengthSquared(const ParticleFlowObject *const pPfo)
{
    if (!LArPfoHelper::IsTwoD(pPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float lengthSquared(0.f);

    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pPfoCluster = *cIter;

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        lengthSquared += LArClusterHelper::GetLengthSquared(pPfoCluster);
    }

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetThreeDLengthSquared(const ParticleFlowObject *const pPfo)
{
    if (!LArPfoHelper::IsThreeD(pPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float lengthSquared(0.f);

    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pPfoCluster = *cIter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        lengthSquared += LArClusterHelper::GetLengthSquared(pPfoCluster);
    }

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetClosestDistance(const ParticleFlowObject *const pPfo, const Cluster *const pCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float bestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pPfoCluster = *iter;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster, pPfoCluster));

        if (thisDistance < bestDistance)
            bestDistance = thisDistance;
    }

    return bestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetTwoDSeparation(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2)
{
    ClusterList clusterListU1, clusterListV1, clusterListW1;
    ClusterList clusterListU2, clusterListV2, clusterListW2;

    LArPfoHelper::GetClusters(pPfo1, TPC_VIEW_U, clusterListU1);
    LArPfoHelper::GetClusters(pPfo1, TPC_VIEW_V, clusterListV1);
    LArPfoHelper::GetClusters(pPfo1, TPC_VIEW_W, clusterListW1);

    LArPfoHelper::GetClusters(pPfo2, TPC_VIEW_U, clusterListU2);
    LArPfoHelper::GetClusters(pPfo2, TPC_VIEW_V, clusterListV2);
    LArPfoHelper::GetClusters(pPfo2, TPC_VIEW_W, clusterListW2);

    float numViews(0.f);
    float distanceSquared(0.f);

    if (!clusterListU1.empty() && !clusterListU2.empty())
    {
        distanceSquared += LArClusterHelper::GetClosestDistance(clusterListU1, clusterListU2);
        numViews += 1.f;
    }

    if (!clusterListV1.empty() && !clusterListV2.empty())
    {
        distanceSquared += LArClusterHelper::GetClosestDistance(clusterListV1, clusterListV2);
        numViews += 1.f;
    }

    if (!clusterListW1.empty() && !clusterListW2.empty())
    {
        distanceSquared += LArClusterHelper::GetClosestDistance(clusterListW1, clusterListW2);
        numViews += 1.f;
    }

    if (numViews < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return std::sqrt(distanceSquared / numViews);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetThreeDSeparation(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2)
{
    ClusterList clusterList1, clusterList2;

    LArPfoHelper::GetClusters(pPfo1, TPC_3D, clusterList1);
    LArPfoHelper::GetClusters(pPfo2, TPC_3D, clusterList2);

    if (clusterList1.empty() || clusterList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return LArClusterHelper::GetClosestDistance(clusterList1, clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsTwoD(const ParticleFlowObject *const pPfo)
{
    for (ClusterList::const_iterator iter = pPfo->GetClusterList().begin(), iterEnd = pPfo->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsThreeD(const ParticleFlowObject *const pPfo)
{
    for (ClusterList::const_iterator iter = pPfo->GetClusterList().begin(), iterEnd = pPfo->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsTrack(const ParticleFlowObject *const pPfo)
{
    const int pdg(pPfo->GetParticleId());

    // muon, pion, proton, kaon
    return ((MU_MINUS == std::abs(pdg)) || (PI_PLUS == std::abs(pdg)) || (PROTON == std::abs(pdg)) || (K_PLUS == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsShower(const ParticleFlowObject *const pPfo)
{
    const int pdg(pPfo->GetParticleId());

    // electron, photon
    return ((E_MINUS == std::abs(pdg)) || (PHOTON == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArPfoHelper::GetPrimaryNeutrino(const ParticleFlowObject *const pPfo)
{
    try
    {
        const ParticleFlowObject *const pParentPfo = LArPfoHelper::GetParentNeutrino(pPfo);
        return pParentPfo->GetParticleId();
    }
    catch (const StatusCodeException &)
    {
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsFinalState(const ParticleFlowObject *const pPfo)
{
    if (pPfo->GetParentPfoList().empty() && !LArPfoHelper::IsNeutrino(pPfo))
        return true;

    if (LArPfoHelper::IsNeutrinoFinalState(pPfo))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsNeutrinoFinalState(const ParticleFlowObject *const pPfo)
{
    return ((pPfo->GetParentPfoList().size() == 1) && (LArPfoHelper::IsNeutrino(*(pPfo->GetParentPfoList().begin()))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::IsNeutrino(const ParticleFlowObject *const pPfo)
{
    const int absoluteParticleId(std::abs(pPfo->GetParticleId()));

    if ((NU_E == absoluteParticleId) || (NU_MU == absoluteParticleId) || (NU_TAU == absoluteParticleId))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetRecoNeutrinos(const PfoList *const pPfoList, PfoList &recoNeutrinos)
{
    if (!pPfoList)
        return;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pPfo(*iter);

        if (LArPfoHelper::IsNeutrino(pPfo))
            recoNeutrinos.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *LArPfoHelper::GetParentPfo(const ParticleFlowObject *const pPfo)
{
    const ParticleFlowObject *pParentPfo = pPfo;

    while (pParentPfo->GetParentPfoList().empty() == false)
    {
        pParentPfo = *(pParentPfo->GetParentPfoList().begin());
    }

    return pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *LArPfoHelper::GetParentNeutrino(const ParticleFlowObject *const pPfo)
{
    const ParticleFlowObject *const pParentPfo = LArPfoHelper::GetParentPfo(pPfo);

    if(!LArPfoHelper::IsNeutrino(pParentPfo))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *LArPfoHelper::GetVertex(const ParticleFlowObject *const pPfo)
{
    if (pPfo->GetVertexList().empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (pPfo->GetVertexList().size() != 1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Vertex *const pVertex = *(pPfo->GetVertexList().begin());

    if (VERTEX_3D != pVertex->GetVertexType())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return pVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetSlidingFitTrajectory(const CartesianPointVector &pointVector, const CartesianVector &vertexPosition,
    const unsigned int layerWindow, const float layerPitch, LArTrackStateVector &trackStateVector)
{
    // Calculate trajectory using 3D sliding fit
    LArTrackTrajectory trackTrajectory;

    try
    {
        const ThreeDSlidingFitResult slidingFitResult(&pointVector, layerWindow, layerPitch);
        const CartesianVector minPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector maxPosition(slidingFitResult.GetGlobalMaxLayerPosition());

        if ((maxPosition - minPosition).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const CartesianVector seedPosition((maxPosition + minPosition) * 0.5f);
        const CartesianVector seedDirection((maxPosition - minPosition).GetUnitVector());

        const float scaleFactor((seedDirection.GetDotProduct(seedPosition - vertexPosition) > 0.f) ? +1.f : -1.f);

        for (CartesianPointVector::const_iterator pIter = pointVector.begin(), pIterEnd =  pointVector.end(); pIter != pIterEnd; ++pIter)
        {
            const CartesianVector nextPoint = *pIter;

            try
            {
                const float rL(slidingFitResult.GetLongitudinalDisplacement(nextPoint));

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
                    LArTrackState(position, direction * scaleFactor, nullptr)));
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

    // Sort trajectory points by distance along track
    std::sort(trackTrajectory.begin(), trackTrajectory.end(), LArPfoHelper::SortByHitProjection);

    for (LArTrackTrajectory::const_iterator tIter = trackTrajectory.begin(), tIterEnd = trackTrajectory.end(); tIter != tIterEnd; ++tIter)
    {
        trackStateVector.push_back(tIter->second);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetSlidingFitTrajectory(const ParticleFlowObject *const pPfo, const Vertex *const pVertex,
    const unsigned int layerWindow, const float layerPitch, LArTrackStateVector &trackStateVector)
{
    // Get 3D vertex position
    const CartesianVector vertexPosition(pVertex->GetPosition());

    // Get 3D hits
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList);

    // Get 3D space points
    CartesianPointVector pointVector;
    for (CaloHitList::const_iterator hIter = caloHitList.begin(), hIterEnd = caloHitList.end(); hIter != hIterEnd; ++hIter)
    {
        const CaloHit *const pCaloHit3D = *hIter;
        pointVector.push_back(pCaloHit3D->GetPositionVector());
    }

    // Check that input list contains some points, and then sort point coordinates by position
    if (pointVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    std::sort(pointVector.begin(), pointVector.end(), LArClusterHelper::SortCoordinatesByPosition);

    // Calculate trajectory using 3D sliding fit
    LArTrackTrajectory trackTrajectory;

    try
    {
        const ThreeDSlidingFitResult slidingFitResult(&pointVector, layerWindow, layerPitch);
        const CartesianVector minPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector maxPosition(slidingFitResult.GetGlobalMaxLayerPosition());

        if ((maxPosition - minPosition).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const CartesianVector seedPosition((maxPosition + minPosition) * 0.5f);
        const CartesianVector seedDirection((maxPosition - minPosition).GetUnitVector());

        const float scaleFactor((seedDirection.GetDotProduct(seedPosition - vertexPosition) > 0.f) ? +1.f : -1.f);

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

    // Sort trajectory points by distance along track
    std::sort(trackTrajectory.begin(), trackTrajectory.end(), LArPfoHelper::SortByHitProjection);

    for (LArTrackTrajectory::const_iterator tIter = trackTrajectory.begin(), tIterEnd = trackTrajectory.end(); tIter != tIterEnd; ++tIter)
    {
        trackStateVector.push_back(tIter->second);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::SortByHitProjection(const LArTrackTrajectoryPoint &lhs, const LArTrackTrajectoryPoint &rhs)
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

bool LArPfoHelper::SortByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    unsigned int nTwoDHitsLhs(0), nThreeDHitsLhs(0); float energyLhs(0.f);
    for (ClusterList::const_iterator iter = pLhs->GetClusterList().begin(), iterEnd = pLhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pClusterLhs = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pClusterLhs))
            nTwoDHitsLhs += pClusterLhs->GetNCaloHits();
        else
            nThreeDHitsLhs += pClusterLhs->GetNCaloHits();

        energyLhs += pClusterLhs->GetHadronicEnergy();
    }

    unsigned int nTwoDHitsRhs(0), nThreeDHitsRhs(0); float energyRhs(0.f);
    for (ClusterList::const_iterator iter = pRhs->GetClusterList().begin(), iterEnd = pRhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pClusterRhs = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pClusterRhs))
            nTwoDHitsRhs += pClusterRhs->GetNCaloHits();
        else
            nThreeDHitsRhs += pClusterRhs->GetNCaloHits();

        energyRhs += pClusterRhs->GetHadronicEnergy();
    }

    if (nTwoDHitsLhs != nTwoDHitsRhs)
        return (nTwoDHitsLhs > nTwoDHitsRhs);

    if (nThreeDHitsLhs != nThreeDHitsRhs)
        return (nThreeDHitsLhs > nThreeDHitsRhs);

    // ATTN Need an efficient (balance with well-motivated) tie-breaker here. Pfo length, for instance, is extremely slow.
    return (energyLhs > energyRhs);
}

} // namespace lar_content
