/**
 *  @file   LArContent/src/LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray shower matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCalculators/LArPseudoLayerCalculator.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayShowerMatchingAlgorithm::Run()
{

    // TODO: Put into tensor format

    this->MergeClusters();


    // TODO: Split clusters

    this->MatchClusters();


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MergeClusters() const
{
    PfoVector inputPfos, seedPfos;
    this->GetPfos(m_inputPfoListName, inputPfos);
    this->SelectPfos(inputPfos, seedPfos);
    this->DeletePfos(seedPfos);

    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    if (clustersU.empty() || clustersV.empty() || clustersW.empty())
        return;

    ParticleList particleList;
    this->MatchViews(clustersU, clustersV, clustersW, particleList);
    this->MergeClusters(particleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MatchClusters() const
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    if (clustersU.empty() || clustersV.empty() || clustersW.empty())
        return;

    ParticleList particleList;
    this->MatchViews(clustersU, clustersV, clustersW, particleList);
    this->MatchClusters(particleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::GetPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
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

void CosmicRayShowerMatchingAlgorithm::GetClusters(const std::string inputClusterListName, ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputClusterListName, pClusterList))

    if (NULL == pClusterList)
    {
        std::cout << "CosmicRayShowerMatchingAlgorithm: could not find cluster list " << inputClusterListName << std::endl;
        return;
    }

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::SelectPfos(const PfoVector &inputPfos, PfoVector &outputPfos) const
{
    for (PfoVector::const_iterator pIter1 = inputPfos.begin(), pIterEnd1 = inputPfos.end(); pIter1 != pIterEnd1; ++pIter1)
    {
        const ParticleFlowObject *pPfo1 = *pIter1;
        const float lengthSquared1(LArPfoHelper::GetTwoDLengthSquared(pPfo1));

        bool isSeed(false);

        for (PfoVector::const_iterator pIter2 = inputPfos.begin(), pIterEnd2 = inputPfos.end(); pIter2 != pIterEnd2; ++pIter2)
        {
            const ParticleFlowObject *pPfo2 = *pIter2;
            const float lengthSquared2(LArPfoHelper::GetTwoDLengthSquared(pPfo2));

            if (pPfo1 == pPfo2)
                continue;

            if (lengthSquared2 < lengthSquared1)
                continue;

            if (LArPfoHelper::GetTwoDSeparation(pPfo1, pPfo2) < m_distanceForMatching)
            {
                isSeed = true;
                break;
            }
        }

        if (isSeed)
            outputPfos.push_back(const_cast<ParticleFlowObject*>(pPfo1));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::DeletePfos(const PfoVector &pfoVector) const
{
    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        ParticleFlowObject *pPfo = *pIter;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, m_inputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MatchViews(const ClusterVector &clusterVectorU, const ClusterVector &clusterVectorV,
    const ClusterVector &clusterVectorW, ParticleList &particleList) const
{
    for (ClusterVector::const_iterator cIterU = clusterVectorU.begin(), cIterEndU = clusterVectorU.end(); cIterU != cIterEndU; ++cIterU)
    {
        const Cluster *pClusterU = *cIterU;

        for (ClusterVector::const_iterator cIterV = clusterVectorV.begin(), cIterEndV = clusterVectorV.end(); cIterV != cIterEndV; ++cIterV)
        {
            const Cluster *pClusterV = *cIterV;

            for (ClusterVector::const_iterator cIterW = clusterVectorW.begin(), cIterEndW = clusterVectorW.end(); cIterW != cIterEndW; ++cIterW)
            {
                const Cluster *pClusterW = *cIterW;

                if (!this->MatchViews(pClusterU, pClusterV, pClusterW))
                    continue;

                particleList.push_back(Particle(pClusterU, pClusterV, pClusterW));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MergeClusters(const ParticleList &particleList) const
{
    ClusterAssociationMap clusterAssociationMap;
    this->IdentifyMerges(particleList, clusterAssociationMap);
    this->ApplyMerges(clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MatchClusters(const ParticleList &inputParticleList) const
{
    ParticleList outputParticleList;
    this->IdentifyMatches(inputParticleList, outputParticleList);
    this->ApplyMatches(outputParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayShowerMatchingAlgorithm::MatchViews(const Cluster *const pClusterU, const Cluster *const pClusterV,
    const Cluster *const pClusterW) const
{
    // Check X overlap
    float xMinU(0.f), xMinV(0.f), xMinW(0.f);
    float xMaxU(0.f), xMaxV(0.f), xMaxW(0.f);
    LArClusterHelper::GetClusterSpanX(pClusterU, xMinU, xMaxU);
    LArClusterHelper::GetClusterSpanX(pClusterV, xMinV, xMaxV);
    LArClusterHelper::GetClusterSpanX(pClusterW, xMinW, xMaxW);

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMinU, std::max(xMinV, xMinW)) - xPitch);
    const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)) + xPitch);
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return false;

    // Match views
    const unsigned int nSamplingPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    for (unsigned int n = 0; n<nSamplingPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nSamplingPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMinU(0.f), zMinV(0.f), zMinW(0.f), zMaxU(0.f), zMaxV(0.f), zMaxW(0.f);
            LArClusterHelper::GetClusterSpanZ(pClusterU, xmin, xmax, zMinU, zMaxU);
            LArClusterHelper::GetClusterSpanZ(pClusterV, xmin, xmax, zMinV, zMaxV);
            LArClusterHelper::GetClusterSpanZ(pClusterW, xmin, xmax, zMinW, zMaxW);

            const float zU(0.5f * (zMinU + zMaxU));
            const float zV(0.5f * (zMinV + zMaxV));
            const float zW(0.5f * (zMinW + zMaxW));

            const float dzU(zMaxU - zMinU);
            const float dzV(zMaxV - zMinV);
            const float dzW(zMaxW - zMinW);
            const float dz(LArGeometryHelper::GetLArPseudoLayerCalculator()->GetZPitch());

            const float zprojU(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, zV, zW));
            const float zprojV(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_W, TPC_VIEW_U, zW, zU));
            const float zprojW(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, zU, zV));

            const float deltaSquared(((zU - zprojU) * (zU - zprojU) + (zV - zprojV) * (zV - zprojV) + (zW - zprojW) * (zW - zprojW)) / 3.f);
            const float sigmaSquared(dzU * dzU + dzV * dzV + dzW * dzW + dz * dz);
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

void CosmicRayShowerMatchingAlgorithm::IdentifyMerges(const ParticleList &particleList, ClusterAssociationMap &clusterAssociationMap) const
{
    for (ParticleList::const_iterator iter1 = particleList.begin(), iterEnd1 = particleList.end(); iter1 != iterEnd1; ++iter1)
    {
        const Particle &particle1 = *iter1;

        for (ParticleList::const_iterator iter2 = iter1, iterEnd2 = particleList.end(); iter2 != iterEnd2; ++iter2)
        {
            const Particle &particle2 = *iter2;

            const bool commonU(particle1.m_pClusterU == particle2.m_pClusterU);
            const bool commonV(particle1.m_pClusterV == particle2.m_pClusterV);
            const bool commonW(particle1.m_pClusterW == particle2.m_pClusterW);

            if (commonU && commonV && commonW)
                continue;

            if (commonU && commonV && !commonW)
                 this->CheckAssociation(particle1.m_pClusterW, particle2.m_pClusterW, clusterAssociationMap);

            if (commonV && commonW && !commonU)
                 this->CheckAssociation(particle1.m_pClusterU, particle2.m_pClusterU, clusterAssociationMap);

            if (commonW && commonU && !commonV)
                 this->CheckAssociation(particle1.m_pClusterV, particle2.m_pClusterV, clusterAssociationMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::IdentifyMatches(const ParticleList &inputParticleList, ParticleList &outputParticleList) const
{
    for (ParticleList::const_iterator iter1 = inputParticleList.begin(), iterEnd1 = inputParticleList.end(); iter1 != iterEnd1; ++iter1)
    {
        const Particle &particle1 = *iter1;

        bool isGoodMatch(true);

        for (ParticleList::const_iterator iter2 = inputParticleList.begin(), iterEnd2 = inputParticleList.end(); iter2 != iterEnd2; ++iter2)
        {
            const Particle &particle2 = *iter2;

            const bool commonU(particle1.m_pClusterU == particle2.m_pClusterU);
            const bool commonV(particle1.m_pClusterV == particle2.m_pClusterV);
            const bool commonW(particle1.m_pClusterW == particle2.m_pClusterW);

            if (commonU && commonV && commonW)
                continue;

            const bool checkU(commonU || 2 * particle2.m_pClusterU->GetNCaloHits() > particle1.m_pClusterU->GetNCaloHits());
            const bool checkV(commonV || 2 * particle2.m_pClusterV->GetNCaloHits() > particle1.m_pClusterV->GetNCaloHits());
            const bool checkW(commonW || 2 * particle2.m_pClusterW->GetNCaloHits() > particle1.m_pClusterW->GetNCaloHits());

            if (checkU && checkV && checkW)
            {
                isGoodMatch = false;
                break;
            }
        }

        if (isGoodMatch)
            outputParticleList.push_back(particle1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::CheckAssociation(const Cluster *const pCluster1, const Cluster *const pCluster2,
    ClusterAssociationMap &clusterAssociationMap) const
{
    if (pCluster1 == pCluster2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (LArClusterHelper::GetClusterHitType(pCluster1) != LArClusterHelper::GetClusterHitType(pCluster2))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (LArClusterHelper::GetClosestDistance(pCluster1, pCluster2) < 10.f)
    {
        clusterAssociationMap[pCluster1].insert(const_cast<Cluster*>(pCluster2));
        clusterAssociationMap[pCluster2].insert(const_cast<Cluster*>(pCluster1));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::ApplyMerges(const ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterAssociationMap clusterMergeMap;
    this->CollectAssociatedClusters(clusterAssociationMap, clusterMergeMap);

    for (ClusterAssociationMap::const_iterator iter1 = clusterMergeMap.begin(), iterEnd1 = clusterMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = const_cast<Cluster*>(iter1->first);
        const ClusterList &mergeList = iter1->second;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pSeedCluster));
        const std::string clusterListName((TPC_VIEW_U == hitType) ? m_inputClusterListNameU :
                                          (TPC_VIEW_V == hitType) ? m_inputClusterListNameV : m_inputClusterListNameW);

        for (ClusterList::iterator iter2 = mergeList.begin(), iterEnd2 = mergeList.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pAssociatedCluster = *iter2;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster,
                clusterListName, clusterListName));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::ApplyMatches(const ParticleList &particleList) const
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

        // TODO - correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
        pfoParameters.m_clusterList = clusterList;

        ParticleFlowObject *pPfo(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    if (!pPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_inputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::CollectAssociatedClusters(const ClusterAssociationMap &clusterAssociationMap,
    ClusterAssociationMap &clusterMergeMap) const
{
    ClusterList vetoList;

    for (ClusterAssociationMap::const_iterator iter1 = clusterAssociationMap.begin(), iterEnd1 = clusterAssociationMap.end();
        iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = const_cast<Cluster*>(iter1->first);

        if (vetoList.count(pSeedCluster))
            continue;

        ClusterList mergeList;
        this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, clusterAssociationMap, vetoList, mergeList);

        if (mergeList.empty())
            continue;

        clusterMergeMap[pSeedCluster].insert(mergeList.begin(), mergeList.end());

        vetoList.insert(mergeList.begin(), mergeList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::CollectAssociatedClusters(Cluster *pSeedCluster, Cluster *pCurrentCluster,
    const ClusterAssociationMap &clusterAssociationMap, const ClusterList &clusterVetoList, ClusterList &associatedClusterList) const
{
    ClusterList::const_iterator iter0 = clusterVetoList.find(pCurrentCluster);

    if (iter0 != clusterVetoList.end())
        return;

    ClusterAssociationMap::const_iterator iter1 = clusterAssociationMap.find(pCurrentCluster);

    if (iter1 == clusterAssociationMap.end())
        return;

    for (ClusterList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        Cluster *pAssociatedCluster = *iter2;

        if (pAssociatedCluster == pSeedCluster)
            continue;

        if (!associatedClusterList.insert(pAssociatedCluster).second)
            continue;

        this->CollectAssociatedClusters(pSeedCluster, pAssociatedCluster, clusterAssociationMap, clusterVetoList, associatedClusterList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayShowerMatchingAlgorithm::Particle::Particle(const Cluster *pClusterU, const Cluster *pClusterV, const Cluster *pClusterW) :
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

StatusCode CosmicRayShowerMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    m_minCaloHitsPerCluster = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_xOverlapWindow = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverlapWindow", m_xOverlapWindow));

    m_distanceForMatching = 5.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceForMatching", m_distanceForMatching));

    m_pseudoChi2Cut = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
