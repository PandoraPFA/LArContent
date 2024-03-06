/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTrackMatchingAlgorithm::CosmicRayTrackMatchingAlgorithm() :
    m_clusterMinLength(10.f), m_vtxXOverlap(3.f), m_minXOverlap(3.f), m_minXOverlapFraction(0.8f), m_maxDisplacement(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SelectCleanClusters(const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    ClusterVector clusterVector;

    // Select long clusters
    for (const Cluster *const pCluster : inputVector)
    {
        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    // Remove long delta rays
    for (const Cluster *const pCluster : clusterVector)
    {
        const float lengthSquared(LArClusterHelper::GetLengthSquared(pCluster));
        CartesianVector innerVertex(0.f, 0.f, 0.f), outerVertex(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pCluster, innerVertex, outerVertex);

        bool isDeltaRay(false);

        for (const Cluster *const pClusterCheck : clusterVector)
        {
            if (pCluster == pClusterCheck)
                continue;

            if ((LArClusterHelper::GetLengthSquared(pClusterCheck) > 10.f * lengthSquared) &&
                (LArClusterHelper::GetClosestDistance(innerVertex, pClusterCheck) < m_vtxXOverlap ||
                    LArClusterHelper::GetClosestDistance(outerVertex, pClusterCheck) < m_vtxXOverlap))
            {
                isDeltaRay = true;
                break;
            }
        }

        if (!isDeltaRay)
            outputVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackMatchingAlgorithm::MatchClusters(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    // Use start and end points
    CartesianVector innerVertex1(0.f, 0.f, 0.f), outerVertex1(0.f, 0.f, 0.f);
    CartesianVector innerVertex2(0.f, 0.f, 0.f), outerVertex2(0.f, 0.f, 0.f);

    LArClusterHelper::GetExtremalCoordinates(pCluster1, innerVertex1, outerVertex1);
    LArClusterHelper::GetExtremalCoordinates(pCluster2, innerVertex2, outerVertex2);

    const float dxA(std::fabs(innerVertex2.GetX() - innerVertex1.GetX())); // inner1 <-> inner2
    const float dxB(std::fabs(outerVertex2.GetX() - outerVertex1.GetX())); // outer1 <-> outer2

    const float dxC(std::fabs(outerVertex2.GetX() - innerVertex1.GetX())); // inner1 <-> outer2
    const float dxD(std::fabs(innerVertex2.GetX() - outerVertex1.GetX())); // inner2 <-> outer1

    const float xVertex(std::min(std::max(dxA, dxB), std::max(dxC, dxD)));

    if (xVertex < m_vtxXOverlap)
        return true;

    // use X overlap
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));
    const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));

    if (xSpan < std::numeric_limits<float>::epsilon())
        return false;

    if (xOverlap > m_minXOverlap && xOverlap / xSpan > m_minXOverlapFraction)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackMatchingAlgorithm::CheckMatchedClusters3D(
    const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3) const
{
    // Check that three clusters have a consistent 3D position
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

    if (hitType1 == hitType2 || hitType2 == hitType3 || hitType3 == hitType1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    CartesianVector innerVertex1(0.f, 0.f, 0.f), outerVertex1(0.f, 0.f, 0.f);
    CartesianVector innerVertex2(0.f, 0.f, 0.f), outerVertex2(0.f, 0.f, 0.f);
    CartesianVector innerVertex3(0.f, 0.f, 0.f), outerVertex3(0.f, 0.f, 0.f);

    LArClusterHelper::GetExtremalCoordinates(pCluster1, innerVertex1, outerVertex1);
    LArClusterHelper::GetExtremalCoordinates(pCluster2, innerVertex2, outerVertex2);
    LArClusterHelper::GetExtremalCoordinates(pCluster3, innerVertex3, outerVertex3);

    for (unsigned int n = 0; n < 4; ++n)
    {
        CartesianVector vtx1(1 == n ? outerVertex1 : innerVertex1);
        CartesianVector end1(1 == n ? innerVertex1 : outerVertex1);

        CartesianVector vtx2(2 == n ? outerVertex2 : innerVertex2);
        CartesianVector end2(2 == n ? innerVertex2 : outerVertex2);

        CartesianVector vtx3(3 == n ? outerVertex3 : innerVertex3);
        CartesianVector end3(3 == n ? innerVertex3 : outerVertex3);

        if (std::fabs(vtx1.GetX() - vtx2.GetX()) < std::max(m_vtxXOverlap, std::fabs(vtx1.GetX() - end2.GetX())) &&
            std::fabs(end1.GetX() - end2.GetX()) < std::max(m_vtxXOverlap, std::fabs(end1.GetX() - vtx2.GetX())) &&
            std::fabs(vtx2.GetX() - vtx3.GetX()) < std::max(m_vtxXOverlap, std::fabs(vtx2.GetX() - end3.GetX())) &&
            std::fabs(end2.GetX() - end3.GetX()) < std::max(m_vtxXOverlap, std::fabs(end2.GetX() - vtx3.GetX())) &&
            std::fabs(vtx3.GetX() - vtx1.GetX()) < std::max(m_vtxXOverlap, std::fabs(vtx3.GetX() - end1.GetX())) &&
            std::fabs(end3.GetX() - end1.GetX()) < std::max(m_vtxXOverlap, std::fabs(end3.GetX() - vtx1.GetX())))
        {
            float chi2(0.f);
            CartesianVector projVtx1(0.f, 0.f, 0.f), projEnd1(0.f, 0.f, 0.f);
            CartesianVector projVtx2(0.f, 0.f, 0.f), projEnd2(0.f, 0.f, 0.f);
            CartesianVector projVtx3(0.f, 0.f, 0.f), projEnd3(0.f, 0.f, 0.f);

            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, vtx1, vtx2, projVtx3, chi2);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, end1, end2, projEnd3, chi2);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, vtx2, vtx3, projVtx1, chi2);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, end2, end3, projEnd1, chi2);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType3, hitType1, vtx3, vtx1, projVtx2, chi2);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType3, hitType1, end3, end1, projEnd2, chi2);

            const bool matchedVtx1(LArClusterHelper::GetClosestDistance(projVtx1, pCluster1) < m_maxDisplacement);
            const bool matchedVtx2(LArClusterHelper::GetClosestDistance(projVtx2, pCluster2) < m_maxDisplacement);
            const bool matchedVtx3(LArClusterHelper::GetClosestDistance(projVtx3, pCluster3) < m_maxDisplacement);

            const bool matchedEnd1(LArClusterHelper::GetClosestDistance(projEnd1, pCluster1) < m_maxDisplacement);
            const bool matchedEnd2(LArClusterHelper::GetClosestDistance(projEnd2, pCluster2) < m_maxDisplacement);
            const bool matchedEnd3(LArClusterHelper::GetClosestDistance(projEnd3, pCluster3) < m_maxDisplacement);

            const bool matchedCluster1(matchedVtx1 || matchedEnd1);
            const bool matchedCluster2(matchedVtx2 || matchedEnd2);
            const bool matchedCluster3(matchedVtx3 || matchedEnd3);
            const bool matchedVtx(matchedVtx1 || matchedVtx2 || matchedVtx3);
            const bool matchedEnd(matchedEnd1 || matchedEnd2 || matchedEnd3);

            if (matchedCluster1 && matchedCluster2 && matchedCluster3 && matchedVtx && matchedEnd)
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SetPfoParameters(const Particle &particle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    ClusterList clusterList;
    if (particle.m_pClusterU)
        clusterList.push_back(particle.m_pClusterU);
    if (particle.m_pClusterV)
        clusterList.push_back(particle.m_pClusterV);
    if (particle.m_pClusterW)
        clusterList.push_back(particle.m_pClusterW);

    // TODO Correct these placeholder parameters
    pfoParameters.m_particleId = MU_MINUS; // TRACK
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterMinLength", m_clusterMinLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VtxXOverlap", m_vtxXOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlap", m_minXOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDisplacement", m_maxDisplacement));

    return CosmicRayBaseMatchingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
