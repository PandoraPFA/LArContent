/**
 *  @file   LArContent/src/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray shower matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeltaRayMatchingAlgorithm::DeltaRayMatchingAlgorithm() :
    m_minCaloHitsPerCluster(3),
    m_xOverlapWindow(1.f),
    m_distanceForMatching(5.f),
    m_pseudoChi2Cut(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMatchingAlgorithm::Run()
{
    PfoVector pfoVector;
    this->GetAllPfos(m_parentPfoListName, pfoVector);

    if (pfoVector.empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "DeltaRayMatchingAlgorithm: pfo list " << m_parentPfoListName << " unavailable." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    this->ThreeViewMatching();
    this->TwoViewMatching();
    this->OneViewMatching();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetAllPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
        pfoVector.push_back(*iter);

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetTrackPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
{
    PfoVector inputVector;
    this->GetAllPfos(inputPfoListName, inputVector);

    for (PfoVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;

        if (!LArPfoHelper::IsTrack(pPfo))
            continue;

        pfoVector.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetClusters(const std::string &inputClusterListName, pandora::ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputClusterListName, pClusterList))

    if (NULL == pClusterList)
        return;

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;

        if (!pCluster->IsAvailable() || (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::ThreeViewMatching() const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    ParticleList initialParticleList, finalParticleList;
    this->ThreeViewMatching(clustersU, clustersV, clustersW, initialParticleList);
    this->SelectParticles(initialParticleList, finalParticleList);
    this->CreateParticles(finalParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::TwoViewMatching() const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    ParticleList initialParticleList, finalParticleList;
    this->TwoViewMatching(clustersU, clustersV, initialParticleList);
    this->TwoViewMatching(clustersV, clustersW, initialParticleList);
    this->TwoViewMatching(clustersW, clustersU, initialParticleList);
    this->SelectParticles(initialParticleList, finalParticleList);
    this->CreateParticles(finalParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::OneViewMatching() const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    ParticleList initialParticleList, finalParticleList;
    this->ThreeViewMatching(clustersU, clustersV, clustersW, initialParticleList);
    this->OneViewMatching(clustersU, initialParticleList);
    this->OneViewMatching(clustersV, initialParticleList);
    this->OneViewMatching(clustersW, initialParticleList);
    this->SelectParticles(initialParticleList, finalParticleList);
    this->CreateParticles(finalParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::ThreeViewMatching(const ClusterVector &clusters1, const ClusterVector &clusters2, const ClusterVector &clusters3,
    ParticleList &particleList) const
{
    if (clusters1.empty() || clusters2.empty() || clusters3.empty())
        return;

    for (ClusterVector::const_iterator cIter1 = clusters1.begin(), cIterEnd1 = clusters1.end(); cIter1 != cIterEnd1; ++cIter1)
    {
        Cluster *const pCluster1 = *cIter1;

        for (ClusterVector::const_iterator cIter2 = clusters2.begin(), cIterEnd2 = clusters2.end(); cIter2 != cIterEnd2; ++cIter2)
        {
            Cluster *const pCluster2 = *cIter2;

            for (ClusterVector::const_iterator cIter3 = clusters3.begin(), cIterEnd3 = clusters3.end(); cIter3 != cIterEnd3; ++cIter3)
            {
                Cluster *const pCluster3 = *cIter3;

                if (!pCluster1->IsAvailable() || !pCluster2->IsAvailable() || !pCluster3->IsAvailable())
                    continue;

                if (!this->AreClustersMatched(pCluster1, pCluster2, pCluster3))
                    continue;

                ParticleFlowObject *pBestPfo = NULL;
                this->FindBestParentPfo(pCluster1, pCluster2, pCluster3, pBestPfo);

                // ATTN Need to record all matches when all three views are used
                particleList.push_back(Particle(pCluster1, pCluster2, pCluster3, pBestPfo));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::TwoViewMatching(const ClusterVector &clusters1, const ClusterVector &clusters2, ParticleList &particleList) const
{
    if (clusters1.empty() || clusters2.empty())
        return;

    for (ClusterVector::const_iterator cIter1 = clusters1.begin(), cIterEnd1 = clusters1.end(); cIter1 != cIterEnd1; ++cIter1)
    {
        Cluster *const pCluster1 = *cIter1;

        for (ClusterVector::const_iterator cIter2 = clusters2.begin(), cIterEnd2 = clusters2.end(); cIter2 != cIterEnd2; ++cIter2)
        {
            Cluster *const pCluster2 = *cIter2;

            if (!pCluster1->IsAvailable() || !pCluster2->IsAvailable())
                continue;

            if (!this->AreClustersMatched(pCluster1, pCluster2, NULL))
                continue;

            ParticleFlowObject *pBestPfo = NULL;
            this->FindBestParentPfo(pCluster1, pCluster2, NULL, pBestPfo);

            if (NULL == pBestPfo)
                continue;

            particleList.push_back(Particle(pCluster1, pCluster2, NULL, pBestPfo));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::OneViewMatching(const ClusterVector &clusters, ParticleList &particleList) const
{
    if (clusters.empty())
        return;

    for (ClusterVector::const_iterator cIter = clusters.begin(), cIterEnd = clusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *const pCluster = *cIter;

        if (!pCluster->IsAvailable())
            continue;

        ParticleFlowObject *pBestPfo = NULL;

        this->FindBestParentPfo(pCluster, NULL, NULL, pBestPfo);

        if (NULL == pBestPfo)
            continue;

        particleList.push_back(Particle(pCluster, NULL, NULL, pBestPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::SelectParticles(const ParticleList &initialParticles, ParticleList &finalParticles) const
{
    for (ParticleList::const_iterator iter1 = initialParticles.begin(), iterEnd1 = initialParticles.end(); iter1 != iterEnd1; ++iter1)
    {
        const Particle &particle1 = *iter1;

        bool isGoodParticle(true);

        for (ParticleList::const_iterator iter2 = initialParticles.begin(), iterEnd2 = initialParticles.end(); iter2 != iterEnd2; ++iter2)
        {
            const Particle &particle2 = *iter2;

            const bool commonU(particle1.GetClusterU() == particle2.GetClusterU());
            const bool commonV(particle1.GetClusterV() == particle2.GetClusterV());
            const bool commonW(particle1.GetClusterW() == particle2.GetClusterW());

            const bool ambiguousU(commonU && NULL != particle1.GetClusterU());
            const bool ambiguousV(commonV && NULL != particle1.GetClusterV());
            const bool ambiguousW(commonW && NULL != particle1.GetClusterW());

            if (commonU && commonV && commonW)
                continue;

            if (ambiguousU || ambiguousV || ambiguousW)
            {
                if (particle2.GetNViews() > particle1.GetNViews())
                {
                    isGoodParticle = false;
                }
                else if (particle2.GetNViews() == particle1.GetNViews() && NULL != particle2.GetParentPfo())
                {
                    if ((NULL == particle1.GetParentPfo()) || (particle2.GetNCaloHits() > particle1.GetNCaloHits()) ||
                        (particle2.GetNCaloHits() == particle1.GetNCaloHits() && particle2.GetLengthSquared() >= particle1.GetLengthSquared()))
                            isGoodParticle = false;
                }

                if (!isGoodParticle)
                    break;
            }
        }

        if (isGoodParticle && NULL != particle1.GetParentPfo())
            finalParticles.push_back(particle1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::CreateParticles(const ParticleList &particleList) const
{
    PfoVector parentVector, daughterVector;
    this->GetTrackPfos(m_parentPfoListName, parentVector);
    this->GetAllPfos(m_daughterPfoListName, daughterVector);

    PfoList parentList(parentVector.begin(), parentVector.end());
    PfoList daughterList(daughterVector.begin(), daughterVector.end());

    for (ParticleList::const_iterator iter = particleList.begin(), iterEnd = particleList.end(); iter != iterEnd; ++iter)
    {
        const Particle &particle = *iter;

        ParticleFlowObject *pParentPfo = particle.GetParentPfo();

        if (NULL == pParentPfo)
            continue;

        Cluster* pClusterU = particle.GetClusterU();
        Cluster* pClusterV = particle.GetClusterV();
        Cluster* pClusterW = particle.GetClusterW();

        if (NULL == pClusterU && NULL == pClusterV && NULL == pClusterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterList clusterList;

        if (NULL != pClusterU)
            clusterList.insert(pClusterU);

        if (NULL != pClusterV)
            clusterList.insert(pClusterV);

        if (NULL != pClusterW)
            clusterList.insert(pClusterW);

        if (parentList.count(pParentPfo))
        {
            this->CreateDaughterPfo(clusterList, pParentPfo);
        }
        else if(daughterList.count(pParentPfo))
        {
            this->AddToDaughterPfo(clusterList, pParentPfo);
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMatchingAlgorithm::AreClustersMatched(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const Cluster *const pCluster3) const
{
    if (NULL == pCluster1 && NULL == pCluster2 && NULL == pCluster3)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // First step: Check X overlap
    float xMin1(-std::numeric_limits<float>::max()), xMax1(+std::numeric_limits<float>::max());
    float xMin2(-std::numeric_limits<float>::max()), xMax2(+std::numeric_limits<float>::max());
    float xMin3(-std::numeric_limits<float>::max()), xMax3(+std::numeric_limits<float>::max());

    if (NULL != pCluster1)
        LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);

    if (NULL != pCluster2)
        LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    if (NULL != pCluster3)
        LArClusterHelper::GetClusterSpanX(pCluster3, xMin3, xMax3);

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, std::max(xMin2, xMin3)) - xPitch);
    const float xMax(std::min(xMax1, std::min(xMax2, xMax3)) + xPitch);
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return false;

    if (NULL == pCluster1 || NULL == pCluster2 || NULL == pCluster3)
        return true;

    // Second step: Check 3D matching
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

    if (hitType1 == hitType2 ||  hitType2 == hitType3 || hitType3 == hitType1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const unsigned int nSamplingPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    for (unsigned int n = 0; n<nSamplingPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nSamplingPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMin1(0.f), zMin2(0.f), zMin3(0.f), zMax1(0.f), zMax2(0.f), zMax3(0.f);
            LArClusterHelper::GetClusterSpanZ(pCluster1, xmin, xmax, zMin1, zMax1);
            LArClusterHelper::GetClusterSpanZ(pCluster2, xmin, xmax, zMin2, zMax2);
            LArClusterHelper::GetClusterSpanZ(pCluster3, xmin, xmax, zMin3, zMax3);

            const float z1(0.5f * (zMin1 + zMax1));
            const float z2(0.5f * (zMin2 + zMax2));
            const float z3(0.5f * (zMin3 + zMax3));

            const float dz1(zMax1 - zMin1);
            const float dz2(zMax2 - zMin2);
            const float dz3(zMax3 - zMin3);
            const float dz4(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

            const float zproj1(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, z2, z3));
            const float zproj2(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType3, hitType1, z3, z1));
            const float zproj3(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, z1, z2));

            const float deltaSquared(((z1 - zproj1) * (z1 - zproj1) + (z2 - zproj2) * (z2 - zproj2) + (z3 - zproj3) * (z3 - zproj3)) / 3.f);
            const float sigmaSquared(dz1 * dz1 + dz2 * dz2 + dz3 * dz3 + dz4 * dz4);
            const float pseudoChi2(deltaSquared / sigmaSquared);

            if (pseudoChi2 < m_pseudoChi2Cut)
                return true;
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::FindBestParentPfo(Cluster *const pCluster1, Cluster *const pCluster2, Cluster *const pCluster3,
    ParticleFlowObject* &pBestPfo) const

{
    PfoVector pfoVector;
    this->GetTrackPfos(m_parentPfoListName, pfoVector);
    this->GetAllPfos(m_daughterPfoListName, pfoVector);

    if (pfoVector.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    float numViews(0.f);
    float lengthSquared(0.f);

    if (NULL != pCluster1)
    {
        lengthSquared += LArClusterHelper::GetLengthSquared(pCluster1);
        numViews += 1.f;
    }

    if (NULL != pCluster2)
    {
        lengthSquared += LArClusterHelper::GetLengthSquared(pCluster2);
        numViews += 1.f;
    }

    if (NULL != pCluster3)
    {
        lengthSquared += LArClusterHelper::GetLengthSquared(pCluster3);
        numViews += 1.f;
    }

    float bestDistanceSquared(numViews * m_distanceForMatching * m_distanceForMatching);

    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;

        if (lengthSquared > LArPfoHelper::GetTwoDLengthSquared(pPfo))
            continue;

        try
        {
            float distanceSquared(0.f);

            if (NULL != pCluster1)
                distanceSquared += this->GetDistanceSquaredToPfo(pCluster1, pPfo);

            if (NULL != pCluster2)
                distanceSquared += this->GetDistanceSquaredToPfo(pCluster2, pPfo);

            if (NULL != pCluster3)
                distanceSquared += this->GetDistanceSquaredToPfo(pCluster3, pPfo);

            if (distanceSquared < bestDistanceSquared)
            {
                pBestPfo = pPfo;
                bestDistanceSquared = distanceSquared;
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode()))
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::GetDistanceSquaredToPfo(const Cluster *const pCluster, const ParticleFlowObject *const pPfo) const
{
    ClusterList pfoClusterList;
    LArPfoHelper::GetClusters(pPfo, LArClusterHelper::GetClusterHitType(pCluster), pfoClusterList);

    return LArClusterHelper::GetClosestDistance(pCluster, pfoClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::CreateDaughterPfo(const ClusterList &clusterList, ParticleFlowObject *const pParentPfo) const
{
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    // TODO Correct these placeholder parameters
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = E_MINUS; // SHOWER
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusterList;

    ParticleFlowObject *pDaughterPfo(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pDaughterPfo));

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_daughterPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::AddToDaughterPfo(const ClusterList &clusterList, ParticleFlowObject *const pParentPfo) const
{
    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pDaughterCluster = *cIter;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const std::string clusterListName((TPC_VIEW_U == hitType) ? m_inputClusterListNameU :
                                          (TPC_VIEW_V == hitType) ? m_inputClusterListNameV : m_inputClusterListNameW);

        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pParentPfo, hitType, pfoClusters);

        if (pfoClusters.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Cluster *pParentCluster = *(pfoClusters.begin());

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
            clusterListName, clusterListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRayMatchingAlgorithm::Particle::Particle(Cluster *const pCluster1, Cluster *const pCluster2, Cluster *const pCluster3,
    ParticleFlowObject *const pPfo) :
    m_pClusterU(NULL),
    m_pClusterV(NULL),
    m_pClusterW(NULL),
    m_pParentPfo(NULL)
{
    const HitType hitType1(NULL != pCluster1 ? LArClusterHelper::GetClusterHitType(pCluster1) : HIT_CUSTOM);
    const HitType hitType2(NULL != pCluster2 ? LArClusterHelper::GetClusterHitType(pCluster2) : HIT_CUSTOM);
    const HitType hitType3(NULL != pCluster3 ? LArClusterHelper::GetClusterHitType(pCluster3) : HIT_CUSTOM);

    m_pClusterU = ((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : (TPC_VIEW_U == hitType3) ? pCluster3 : NULL);
    m_pClusterV = ((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : (TPC_VIEW_V == hitType3) ? pCluster3 : NULL);
    m_pClusterW = ((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : (TPC_VIEW_W == hitType3) ? pCluster3 : NULL);
    m_pParentPfo = pPfo;

    if (NULL == m_pClusterU && NULL == m_pClusterV && NULL == m_pClusterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DeltaRayMatchingAlgorithm::Particle::GetNViews() const
{
    unsigned int numViews(0);

    if (NULL != m_pClusterU)
        numViews += 1;

    if (NULL != m_pClusterV)
        numViews += 1;

    if (NULL != m_pClusterW)
        numViews += 1;

    return numViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DeltaRayMatchingAlgorithm::Particle::GetNCaloHits() const
{
    unsigned int numCaloHits(0);

    if (NULL != m_pClusterU)
        numCaloHits += m_pClusterU->GetNCaloHits();

    if (NULL != m_pClusterV)
        numCaloHits += m_pClusterV->GetNCaloHits();

    if (NULL != m_pClusterW)
        numCaloHits += m_pClusterW->GetNCaloHits();

    return numCaloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::Particle::GetLengthSquared() const
{
    float lengthSquared(0.f);

    if (NULL != m_pClusterU)
        lengthSquared += LArClusterHelper::GetLengthSquared(m_pClusterU);

    if (NULL != m_pClusterV)
        lengthSquared += LArClusterHelper::GetLengthSquared(m_pClusterV);

    if (NULL != m_pClusterW)
        lengthSquared += LArClusterHelper::GetLengthSquared(m_pClusterW);

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentPfoListName", m_parentPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterPfoListName", m_daughterPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverlapWindow", m_xOverlapWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceForMatching", m_distanceForMatching));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
