/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayBaseMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayBaseMatchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

StatusCode CosmicRayBaseMatchingAlgorithm::Run()
{
    // Get the available clusters for each view
    ClusterVector availableClustersU, availableClustersV, availableClustersW;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNameU, availableClustersU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNameV, availableClustersV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNameW, availableClustersW));

    // Select clean clusters in each view
    ClusterVector cleanClustersU, cleanClustersV, cleanClustersW;
    this->SelectCleanClusters(availableClustersU, cleanClustersU);
    this->SelectCleanClusters(availableClustersV, cleanClustersV);
    this->SelectCleanClusters(availableClustersW, cleanClustersW);

    // Build associations between pairs of views
    ClusterAssociationMap matchedClusterUV, matchedClusterVW, matchedClusterWU;
    this->MatchClusters(cleanClustersU, cleanClustersV, matchedClusterUV);
    this->MatchClusters(cleanClustersV, cleanClustersW, matchedClusterVW);
    this->MatchClusters(cleanClustersW, cleanClustersU, matchedClusterWU);

    // Build particles from associations
    ParticleList particleList;
    this->MatchThreeViews(matchedClusterUV, matchedClusterVW, matchedClusterWU, particleList);
    this->MatchTwoViews(matchedClusterUV, matchedClusterVW, matchedClusterWU, particleList);
    this->BuildParticles(particleList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayBaseMatchingAlgorithm::GetAvailableClusters(const std::string inputClusterListName, ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, inputClusterListName, pClusterList))

    if (!pClusterList || pClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CosmicRayBaseMatchingAlgorithm: unable to find cluster list " << inputClusterListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBaseMatchingAlgorithm::MatchClusters(
    const ClusterVector &clusterVector1, const ClusterVector &clusterVector2, ClusterAssociationMap &matchedClusters12) const
{
    // Check that there are input clusters from both views
    if (clusterVector1.empty() || clusterVector2.empty())
        return;

    const HitType hitType1(LArClusterHelper::GetClusterHitType(*clusterVector1.begin()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(*clusterVector2.begin()));

    if (hitType1 == hitType2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (const Cluster *const pCluster1 : clusterVector1)
    {
        for (const Cluster *const pCluster2 : clusterVector2)
        {
            if (this->MatchClusters(pCluster1, pCluster2))
            {
                UIntSet daughterVolumeIntersection;
                LArGeometryHelper::GetCommonDaughterVolumes(pCluster1, pCluster2, daughterVolumeIntersection);

                if (!daughterVolumeIntersection.empty())
                    matchedClusters12[pCluster1].push_back(pCluster2);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBaseMatchingAlgorithm::MatchThreeViews(const ClusterAssociationMap &matchedClusters12,
    const ClusterAssociationMap &matchedClusters23, const ClusterAssociationMap &matchedClusters31, ParticleList &matchedParticles) const
{
    if (matchedClusters12.empty() || matchedClusters23.empty() || matchedClusters31.empty())
        return;

    ParticleList candidateParticles;

    ClusterList clusterList1;
    for (const auto &mapEntry : matchedClusters12)
        clusterList1.push_back(mapEntry.first);
    clusterList1.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : clusterList1)
    {
        const ClusterList &clusterList2(matchedClusters12.at(pCluster1));

        for (const Cluster *const pCluster2 : clusterList2)
        {
            ClusterAssociationMap::const_iterator iter23 = matchedClusters23.find(pCluster2);

            if (matchedClusters23.end() == iter23)
                continue;

            const ClusterList &clusterList3 = iter23->second;

            for (const Cluster *const pCluster3 : clusterList3)
            {
                ClusterAssociationMap::const_iterator iter31 = matchedClusters31.find(pCluster3);

                if (matchedClusters31.end() == iter31)
                    continue;

                if (iter31->second.end() == std::find(iter31->second.begin(), iter31->second.end(), pCluster1))
                    continue;

                const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
                const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
                const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

                if (!this->CheckMatchedClusters3D(pCluster1, pCluster2, pCluster3))
                    continue;

                const Cluster *const pClusterU((TPC_VIEW_U == hitType1)   ? pCluster1
                                               : (TPC_VIEW_U == hitType2) ? pCluster2
                                               : (TPC_VIEW_U == hitType3) ? pCluster3
                                                                          : NULL);
                const Cluster *const pClusterV((TPC_VIEW_V == hitType1)   ? pCluster1
                                               : (TPC_VIEW_V == hitType2) ? pCluster2
                                               : (TPC_VIEW_V == hitType3) ? pCluster3
                                                                          : NULL);
                const Cluster *const pClusterW((TPC_VIEW_W == hitType1)   ? pCluster1
                                               : (TPC_VIEW_W == hitType2) ? pCluster2
                                               : (TPC_VIEW_W == hitType3) ? pCluster3
                                                                          : NULL);

                candidateParticles.push_back(Particle(pClusterU, pClusterV, pClusterW));
            }
        }
    }

    return this->ResolveAmbiguities(candidateParticles, matchedParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBaseMatchingAlgorithm::MatchTwoViews(const ClusterAssociationMap &matchedClusters12,
    const ClusterAssociationMap &matchedClusters23, const ClusterAssociationMap &matchedClusters31, ParticleList &matchedParticles) const
{
    ParticleList candidateParticles;
    this->MatchTwoViews(matchedClusters12, candidateParticles);
    this->MatchTwoViews(matchedClusters23, candidateParticles);
    this->MatchTwoViews(matchedClusters31, candidateParticles);

    return this->ResolveAmbiguities(candidateParticles, matchedParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBaseMatchingAlgorithm::MatchTwoViews(const ClusterAssociationMap &matchedClusters12, ParticleList &matchedParticles) const
{
    if (matchedClusters12.empty())
        return;

    ClusterList clusterList1;
    for (const auto &mapEntry : matchedClusters12)
        clusterList1.push_back(mapEntry.first);
    clusterList1.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : clusterList1)
    {
        const ClusterList &clusterList2(matchedClusters12.at(pCluster1));

        for (const Cluster *const pCluster2 : clusterList2)
        {
            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            const Cluster *const pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
            const Cluster *const pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
            const Cluster *const pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);

            matchedParticles.push_back(Particle(pClusterU, pClusterV, pClusterW));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBaseMatchingAlgorithm::ResolveAmbiguities(const ParticleList &candidateParticles, ParticleList &matchedParticles) const
{
    for (const Particle &particle1 : candidateParticles)
    {
        bool isGoodMatch(true);

        for (const Particle &particle2 : candidateParticles)
        {
            const bool commonU(particle1.m_pClusterU == particle2.m_pClusterU);
            const bool commonV(particle1.m_pClusterV == particle2.m_pClusterV);
            const bool commonW(particle1.m_pClusterW == particle2.m_pClusterW);

            const bool ambiguousU(commonU && NULL != particle1.m_pClusterU);
            const bool ambiguousV(commonV && NULL != particle1.m_pClusterV);
            const bool ambiguousW(commonW && NULL != particle1.m_pClusterW);

            if (commonU && commonV && commonW)
                continue;

            if (ambiguousU || ambiguousV || ambiguousW)
            {
                isGoodMatch = false;
                break;
            }
        }

        if (isGoodMatch)
            matchedParticles.push_back(particle1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBaseMatchingAlgorithm::BuildParticles(const ParticleList &particleList)
{
    if (particleList.empty())
        return;

    const PfoList *pPfoList = NULL;
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Particle &particle : particleList)
    {
        const Cluster *const pClusterU = particle.m_pClusterU;
        const Cluster *const pClusterV = particle.m_pClusterV;
        const Cluster *const pClusterW = particle.m_pClusterW;

        const bool isAvailableU((NULL != pClusterU) ? pClusterU->IsAvailable() : true);
        const bool isAvailableV((NULL != pClusterV) ? pClusterV->IsAvailable() : true);
        const bool isAvailableW((NULL != pClusterW) ? pClusterW->IsAvailable() : true);

        if (!(isAvailableU && isAvailableV && isAvailableW))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        this->SetPfoParameters(particle, pfoParameters);

        if (pfoParameters.m_clusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const ParticleFlowObject *pPfo(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    if (!pPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayBaseMatchingAlgorithm::Particle::Particle(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW) :
    m_pClusterU(pClusterU), m_pClusterV(pClusterV), m_pClusterW(pClusterW)
{
    if (NULL == m_pClusterU && NULL == m_pClusterV && NULL == m_pClusterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const HitType hitTypeU(NULL == m_pClusterU ? TPC_VIEW_U : LArClusterHelper::GetClusterHitType(m_pClusterU));
    const HitType hitTypeV(NULL == m_pClusterV ? TPC_VIEW_V : LArClusterHelper::GetClusterHitType(m_pClusterV));
    const HitType hitTypeW(NULL == m_pClusterW ? TPC_VIEW_W : LArClusterHelper::GetClusterHitType(m_pClusterW));

    if (!(TPC_VIEW_U == hitTypeU && TPC_VIEW_V == hitTypeV && TPC_VIEW_W == hitTypeW))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayBaseMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
