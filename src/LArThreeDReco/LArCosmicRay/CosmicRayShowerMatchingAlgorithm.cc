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
    

  

    this->MergeClusters();

    this->MatchClusters();

 
    

   

// --- BEGIN EVENT DISPLAY ---  
// ClusterVector newClustersU, newClustersV, newClustersW;
// this->GetClusters(m_inputClusterListNameU, newClustersU);
// this->GetClusters(m_inputClusterListNameV, newClustersV);
// this->GetClusters(m_inputClusterListNameW, newClustersW);
// ClusterList newListU, newListV, newListW;
// newListU.insert(newClustersU.begin(), newClustersU.end());
// newListV.insert(newClustersV.begin(), newClustersV.end());
// newListW.insert(newClustersW.begin(), newClustersW.end());
// PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
// PandoraMonitoringApi::VisualizeClusters(&newListU, "ClustersU", RED);
// PandoraMonitoringApi::VisualizeClusters(&newListV, "ClustersV", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&newListW, "ClustersW", GREEN);
// PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---



    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MergeClusters()
{
    // Get current pfos and clusters
    PfoVector inputPfos, seedPfos;
    ClusterToPfoMap clusterToPfoMap;
    this->GetPfos(m_inputPfoListName, inputPfos);
    this->SelectPfos(inputPfos, seedPfos);
    this->GetPfoClusterMap(seedPfos, clusterToPfoMap);

    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);

    LArPfoHelper::GetClusters(seedPfos, TPC_VIEW_U, clustersU);
    LArPfoHelper::GetClusters(seedPfos, TPC_VIEW_V, clustersV);
    LArPfoHelper::GetClusters(seedPfos, TPC_VIEW_W, clustersW);

    if (clustersU.empty() || clustersV.empty() || clustersW.empty())
        return;

    ParticleList particleList;
    this->MatchViews(seedPfos, particleList);
    this->MatchViews(clustersU, clustersV, clustersW, particleList);

    ClusterAssociationMap clusterAssociationMap;
    this->BuildAssociationMap(particleList, clusterAssociationMap);
    this->MergeClusters(clusterToPfoMap, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::MatchClusters()
{
    ClusterVector clustersU, clustersV, clustersW;
    this->GetClusters(m_inputClusterListNameU, clustersU);
    this->GetClusters(m_inputClusterListNameV, clustersV);
    this->GetClusters(m_inputClusterListNameW, clustersW);    

    if (clustersU.empty() || clustersV.empty() || clustersW.empty())
        return;

    ParticleList particleList;
    this->MatchViews(clustersU, clustersV, clustersW, particleList);
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

void CosmicRayShowerMatchingAlgorithm::GetPfoClusterMap(const PfoVector &pfoVector, ClusterToPfoMap &clusterToPfoMap) const
{
    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;
        const ClusterList &pfoClusterList = pPfo->GetClusterList();
        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pPfoCluster = *cIter;

            if (TPC_3D == LArClusterHelper::GetClusterHitType(pPfoCluster))
                continue;  

            clusterToPfoMap.insert(ClusterToPfoMap::value_type(pPfoCluster, pPfo));
        }
    }
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

            if (LArPfoHelper::GetTwoDSeparation(pPfo1, pPfo2) < 5.f)
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
 
void CosmicRayShowerMatchingAlgorithm::MatchViews(const PfoVector &pfoVector, ParticleList &particleList) const
{
    ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;
    
    LArPfoHelper::GetClusters(pfoVector, TPC_VIEW_U, clusterVectorU);
    LArPfoHelper::GetClusters(pfoVector, TPC_VIEW_V, clusterVectorV);
    LArPfoHelper::GetClusters(pfoVector, TPC_VIEW_W, clusterVectorW); 

    for (ClusterVector::const_iterator cIterU = clusterVectorU.begin(), cIterEndU = clusterVectorU.end(); cIterU != cIterEndU; ++cIterU)
    {
        const Cluster *pClusterU = *cIterU;

        for (ClusterVector::const_iterator cIterV = clusterVectorV.begin(), cIterEndV = clusterVectorV.end(); cIterV != cIterEndV; ++cIterV)
        {
            const Cluster *pClusterV = *cIterV;

            for (ClusterVector::const_iterator cIterW = clusterVectorW.begin(), cIterEndW = clusterVectorW.end(); cIterW != cIterEndW; ++cIterW)
            {
                const Cluster *pClusterW = *cIterW;

                particleList.push_back(Particle(pClusterU, pClusterV, pClusterW));
            }
        }
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

void CosmicRayShowerMatchingAlgorithm::BuildAssociationMap(const ParticleList &particleList, ClusterAssociationMap &clusterAssociationMap) const
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
                 this->BuildAssociationMap(particle1.m_pClusterW, particle2.m_pClusterW, clusterAssociationMap);

            if (commonV && commonW && !commonU)
                 this->BuildAssociationMap(particle1.m_pClusterU, particle2.m_pClusterU, clusterAssociationMap);

            if (commonW && commonU && !commonV)
                 this->BuildAssociationMap(particle1.m_pClusterV, particle2.m_pClusterV, clusterAssociationMap);
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMatchingAlgorithm::BuildAssociationMap(const Cluster *const pCluster1, const Cluster *const pCluster2, 
    ClusterAssociationMap &clusterAssociationMap) const
{
    if (pCluster1 == pCluster2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (LArClusterHelper::GetClusterHitType(pCluster1) != LArClusterHelper::GetClusterHitType(pCluster2))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (LArClusterHelper::GetClosestDistance(pCluster1, pCluster2) < 5.f)
    {
        clusterAssociationMap[pCluster1].insert(const_cast<Cluster*>(pCluster2));
        clusterAssociationMap[pCluster2].insert(const_cast<Cluster*>(pCluster1));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
void CosmicRayShowerMatchingAlgorithm::MergeClusters(const ClusterToPfoMap &clusterToPfoMap, const ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterAssociationMap clusterMergeMap;
    this->CollectAssociatedClusters(clusterAssociationMap, clusterMergeMap);

    ClusterVector clusterVector;
    
    for (ClusterAssociationMap::const_iterator iter = clusterMergeMap.begin(), iterEnd = clusterMergeMap.end(); iter != iterEnd; ++iter)
    {
        Cluster *pSeedCluster = const_cast<Cluster*>(iter->first);
        const ClusterList &mergeList = iter->second;
        clusterVector.push_back(pSeedCluster);
        clusterVector.insert(clusterVector.end(), mergeList.begin(), mergeList.end());
    }

    PfoList pfoList;

    for (ClusterVector::const_iterator iter1 = clusterVector.begin(), iterEnd1 = clusterVector.end(); iter1 != iterEnd1; ++iter1)
    {
        ClusterToPfoMap::const_iterator iter2 = clusterToPfoMap.find(*iter1);
        if (iter2 == clusterToPfoMap.end())
            continue;

        ParticleFlowObject *pPfo = const_cast<ParticleFlowObject*>(iter2->second);
        pfoList.insert(pPfo);
    }
    
    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, m_inputPfoListName));
    }

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

    m_distanceForMatching = 2.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceForMatching", m_distanceForMatching));

    m_pseudoChi2Cut = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
