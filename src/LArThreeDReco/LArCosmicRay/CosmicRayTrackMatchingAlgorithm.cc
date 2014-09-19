/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTrackMatchingAlgorithm::CosmicRayTrackMatchingAlgorithm() :
    m_clusterMinLength(10.f),
    m_vtxXOverlap(3.f),
    m_minXOverlap(3.f),
    m_minXOverlapFraction(0.8f),
    m_maxDisplacement(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackMatchingAlgorithm::Run()
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

StatusCode CosmicRayTrackMatchingAlgorithm::GetAvailableClusters(const std::string inputClusterListName, ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputClusterListName, pClusterList))

    if (NULL == pClusterList)
    {
        std::cout << "CosmicRayTrackMatchingAlgorithm: could not find cluster list " << inputClusterListName << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;
        if (!pCluster->IsAvailable())
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SelectCleanClusters(const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    ClusterVector clusterVector;

    // Select long clusters
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    // Remove long delta rays
    for (ClusterVector::const_iterator iter1 = clusterVector.begin(), iterEnd1 = clusterVector.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pCluster = *iter1;
        const float lengthSquared(LArClusterHelper::GetLengthSquared(pCluster));
        CartesianVector innerVertex(0.f,0.f,0.f), outerVertex(0.f,0.f,0.f);
        LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, innerVertex, outerVertex);

        bool isDeltaRay(false);

        for (ClusterVector::const_iterator iter2 = clusterVector.begin(), iterEnd2 = clusterVector.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pClusterCheck = *iter2;

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

void CosmicRayTrackMatchingAlgorithm::MatchClusters(const ClusterVector &clusterVector1, const ClusterVector &clusterVector2,
    ClusterAssociationMap &matchedClusters12) const
{
    // Check that there are input clusters from both views
    if (clusterVector1.empty() || clusterVector2.empty())
        return;

    const HitType hitType1(LArClusterHelper::GetClusterHitType(*clusterVector1.begin()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(*clusterVector2.begin()));

    if (hitType1 == hitType2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster* pCluster1 = *iter1;

        for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEnd2 = clusterVector2.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster* pCluster2 = *iter2;

            if (this->MatchClusters(pCluster1, pCluster2))
            {
                matchedClusters12[pCluster1].insert(pCluster2);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackMatchingAlgorithm::MatchClusters(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    // Use start and end points
    CartesianVector innerVertex1(0.f,0.f,0.f), outerVertex1(0.f,0.f,0.f);
    CartesianVector innerVertex2(0.f,0.f,0.f), outerVertex2(0.f,0.f,0.f);

    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster1, innerVertex1, outerVertex1);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster2, innerVertex2, outerVertex2);

    const float dxA(std::fabs(innerVertex2.GetX() - innerVertex1.GetX())); // inner1 <-> inner2
    const float dxB(std::fabs(outerVertex2.GetX() - outerVertex1.GetX())); // outer1 <-> outer2

    const float dxC(std::fabs(outerVertex2.GetX() - innerVertex1.GetX())); // inner1 <-> outer2
    const float dxD(std::fabs(innerVertex2.GetX() - outerVertex1.GetX())); // inner2 <-> outer1

    const float xVertex(std::min(std::max(dxA, dxB), std::max(dxC, dxD)));

    if (xVertex < m_vtxXOverlap)
        return true;

    // use X overlap
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1,xMax2) - std::max(xMin1,xMin2));
    const float xSpan(std::max(xMax1,xMax2) - std::min(xMin1,xMin2));

    if (xOverlap > m_minXOverlap && xOverlap/xSpan > m_minXOverlapFraction)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchThreeViews(const ClusterAssociationMap &matchedClusters12,
    const ClusterAssociationMap &matchedClusters23,  const ClusterAssociationMap &matchedClusters31, ParticleList &matchedParticles) const
{
    if (matchedClusters12.empty() || matchedClusters23.empty() || matchedClusters31.empty())
        return;

    ParticleList candidateParticles;

    for (ClusterAssociationMap::const_iterator iter12 = matchedClusters12.begin(), iterEnd12 = matchedClusters12.end(); iter12 != iterEnd12;
        ++iter12)
    {
        const Cluster* pCluster1 = iter12->first;
        const ClusterList &clusterList2 = iter12->second;

        for(ClusterList::const_iterator iter2 = clusterList2.begin(), iterEnd2 = clusterList2.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pCluster2 = *iter2;

            ClusterAssociationMap::const_iterator iter23 = matchedClusters23.find(pCluster2);
            if (matchedClusters23.end() == iter23)
                continue;

            const ClusterList &clusterList3 = iter23->second;

            for(ClusterList::const_iterator iter3 = clusterList3.begin(), iterEnd3 = clusterList3.end(); iter3 != iterEnd3; ++iter3)
            {
                const Cluster* pCluster3 = *iter3;

                ClusterAssociationMap::const_iterator iter31 = matchedClusters31.find(pCluster3);
                if (matchedClusters31.end() == iter31)
                    continue;

                const ClusterList &clusterList1 = iter31->second;
                Cluster* pCluster1check = const_cast<Cluster*>(pCluster1);
                ClusterList::const_iterator iter1 = clusterList1.find(pCluster1check);

                if (clusterList1.end() == iter1)
                    continue;

                const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
                const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
                const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

                if (!this->CheckMatchedClusters3D(pCluster1, pCluster2, pCluster3))
                    continue;

                const Cluster* pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 :
                                         (TPC_VIEW_U == hitType3) ? pCluster3 : NULL);
                const Cluster* pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 :
                                         (TPC_VIEW_V == hitType3) ? pCluster3 : NULL);
                const Cluster* pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 :
                                         (TPC_VIEW_W == hitType3) ? pCluster3 : NULL);

                candidateParticles.push_back(Particle(pClusterU, pClusterV, pClusterW));
            }
        }
    }

    return this->ResolveAmbiguities(candidateParticles, matchedParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchTwoViews(const ClusterAssociationMap &matchedClusters12,
    const ClusterAssociationMap &matchedClusters23, const ClusterAssociationMap &matchedClusters31, ParticleList &matchedParticles) const
{
    ParticleList candidateParticles;
    this->MatchTwoViews(matchedClusters12, candidateParticles);
    this->MatchTwoViews(matchedClusters23, candidateParticles);
    this->MatchTwoViews(matchedClusters31, candidateParticles);

    return this->ResolveAmbiguities(candidateParticles, matchedParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchTwoViews(const ClusterAssociationMap &matchedClusters12, ParticleList &matchedParticles) const
{
    if (matchedClusters12.empty())
        return;

    for (ClusterAssociationMap::const_iterator iter12 = matchedClusters12.begin(), iterEnd12 = matchedClusters12.end(); iter12 != iterEnd12;
        ++iter12)
    {
        const Cluster* pCluster1 = iter12->first;
        const ClusterList &clusterList2 = iter12->second;

        for (ClusterList::const_iterator iter2 = clusterList2.begin(), iterEnd2 = clusterList2.end() ; iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pCluster2 = *iter2;

            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            const Cluster* pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
            const Cluster* pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
            const Cluster* pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);

            matchedParticles.push_back(Particle(pClusterU, pClusterV, pClusterW));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::ResolveAmbiguities(const ParticleList &candidateParticles, ParticleList &matchedParticles) const
{
    for (ParticleList::const_iterator iter1 = candidateParticles.begin(), iterEnd1 = candidateParticles.end(); iter1 != iterEnd1; ++iter1)
    {
        const Particle &particle1 = *iter1;

        bool isGoodMatch(true);

        for (ParticleList::const_iterator iter2 = candidateParticles.begin(), iterEnd2 = candidateParticles.end(); iter2 != iterEnd2; ++iter2)
        {
            const Particle &particle2 = *iter2;

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

bool CosmicRayTrackMatchingAlgorithm::CheckMatchedClusters3D(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const Cluster *const pCluster3) const
{
    // Check that three clusters have a consistent 3D position
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

    if (hitType1 == hitType2 || hitType2 == hitType3 || hitType3 == hitType1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    CartesianVector innerVertex1(0.f,0.f,0.f), outerVertex1(0.f,0.f,0.f);
    CartesianVector innerVertex2(0.f,0.f,0.f), outerVertex2(0.f,0.f,0.f);
    CartesianVector innerVertex3(0.f,0.f,0.f), outerVertex3(0.f,0.f,0.f);

    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster1, innerVertex1, outerVertex1);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster2, innerVertex2, outerVertex2);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster3, innerVertex3, outerVertex3);

    for (unsigned int n=0; n<4; ++n)
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
            CartesianVector projVtx1(0.f,0.f,0.f), projEnd1(0.f,0.f,0.f);
            CartesianVector projVtx2(0.f,0.f,0.f), projEnd2(0.f,0.f,0.f);
            CartesianVector projVtx3(0.f,0.f,0.f), projEnd3(0.f,0.f,0.f);

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

void CosmicRayTrackMatchingAlgorithm::BuildParticles(const ParticleList &particleList)
{
    if (particleList.empty())
        return;

    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (ParticleList::const_iterator iter = particleList.begin(), iterEnd = particleList.end(); iter != iterEnd; ++iter)
    {
        const Particle &particle = *iter;

        ClusterList clusterList;
        Cluster *pClusterU = const_cast<Cluster*>(particle.m_pClusterU);
        Cluster *pClusterV = const_cast<Cluster*>(particle.m_pClusterV);
        Cluster *pClusterW = const_cast<Cluster*>(particle.m_pClusterW);

        const bool isAvailableU((NULL != pClusterU) ? pClusterU->IsAvailable() : true);
        const bool isAvailableV((NULL != pClusterV) ? pClusterV->IsAvailable() : true);
        const bool isAvailableW((NULL != pClusterW) ? pClusterW->IsAvailable() : true);

        if(!(isAvailableU && isAvailableV && isAvailableW))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (pClusterU) clusterList.insert(pClusterU);
        if (pClusterV) clusterList.insert(pClusterV);
        if (pClusterW) clusterList.insert(pClusterW);

        // TODO Correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = MU_MINUS; // TRACK
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_clusterList = clusterList;

        ParticleFlowObject *pPfo(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    if (!pPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayTrackMatchingAlgorithm::Particle::Particle(const Cluster *pClusterU, const Cluster *pClusterV, const Cluster *pClusterW) :
    m_pClusterU(pClusterU),
    m_pClusterV(pClusterV),
    m_pClusterW(pClusterW)
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

StatusCode CosmicRayTrackMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VtxXOverlap", m_vtxXOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDisplacement", m_maxDisplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
