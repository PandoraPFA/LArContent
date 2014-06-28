/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the cosimic ray shower growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{

void CosmicRayGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable() || (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayGrowingAlgorithm::GetPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
        pfoVector.push_back(*pIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayGrowingAlgorithm::GetPfoClusters(const PfoVector &pfoVector, const HitType hitType, ClusterList &clusterList) const
{
    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;
        this->GetPfoClusters(pPfo, hitType, clusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayGrowingAlgorithm::GetPfoClusters(const ParticleFlowObject *pPfo, const HitType hitType, ClusterList &clusterList) const
{
    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pPfoCluster = *cIter;

        if (hitType != LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        clusterList.insert(pPfoCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryPfoListName", m_primaryPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SecondaryPfoListName", m_secondaryPfoListName));

    m_minCaloHitsPerCluster = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
