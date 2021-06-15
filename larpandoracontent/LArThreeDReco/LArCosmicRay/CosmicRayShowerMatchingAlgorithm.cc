/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

CosmicRayShowerMatchingAlgorithm::CosmicRayShowerMatchingAlgorithm() :
    m_minCaloHitsPerCluster(10),
    m_minXOverlap(1.f),
    m_minXOverlapFraction(0.5f),
    m_pseudoChi2Cut(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::SelectCleanClusters(const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    for (const Cluster *const pCluster : inputVector)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        outputVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayShowerMatchingAlgorithm::MatchClusters(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
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

bool CosmicRayShowerMatchingAlgorithm::CheckMatchedClusters3D(
    const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3) const
{
    // Check that three clusters have a consistent 3D position
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

    if (hitType1 == hitType2 || hitType2 == hitType3 || hitType3 == hitType1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // Requirements on X matching
    float xMin1(0.f), xMin2(0.f), xMin3(0.f), xMax1(0.f), xMax2(0.f), xMax3(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);
    pCluster3->GetClusterSpanX(xMin3, xMax3);

    const float xMin(std::max(xMin1, std::max(xMin2, xMin3)));
    const float xMax(std::min(xMax1, std::min(xMax2, xMax3)));
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return false;

    float p1(std::numeric_limits<float>::max()), p2(std::numeric_limits<float>::max()), p3(std::numeric_limits<float>::max());

    if ((STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pCluster1, xMin, xMax, p1)) ||
        (STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pCluster2, xMin, xMax, p2)) ||
        (STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pCluster3, xMin, xMax, p3)))
    {
        return false;
    }

    const float q3(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, p1, p2));
    const float q1(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, p2, p3));
    const float q2(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType3, hitType1, p3, p1));

    const float pseudoChi2(((q1 - p1) * (q1 - p1) + (q2 - p2) * (q2 - p2) + (q3 - p3) * (q3 - p3)) / 3.f);

    if (pseudoChi2 > m_pseudoChi2Cut)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::SetPfoParameters(const Particle &particle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    ClusterList clusterList;
    if (particle.m_pClusterU)
        clusterList.push_back(particle.m_pClusterU);
    if (particle.m_pClusterV)
        clusterList.push_back(particle.m_pClusterV);
    if (particle.m_pClusterW)
        clusterList.push_back(particle.m_pClusterW);

    // TODO Correct these placeholder parameters
    pfoParameters.m_particleId = E_MINUS; // SHOWER
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlap", m_minXOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PseudoChi2Cut", m_pseudoChi2Cut));

    return CosmicRayBaseMatchingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
