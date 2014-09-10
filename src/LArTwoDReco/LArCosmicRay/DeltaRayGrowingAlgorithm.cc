/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void DeltaRayGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
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

void DeltaRayGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    if (inputClusters.empty())
        return;

    // Get hit type
    const Cluster* pFirstCluster = *(inputClusters.begin());
    const HitType clusterHitType(LArClusterHelper::GetClusterHitType(pFirstCluster));

    // Get parent and daughter Pfos
    PfoVector parentPfos, daughterPfos;
    this->GetPfos(m_parentPfoListName, parentPfos);
    this->GetPfos(m_daughterPfoListName, daughterPfos);

    // TODO Think about sorting of the seed cluster list; currently pfos with most hits contribute first, but their daughter clusters are unsorted
    ClusterList parentClusters, daughterClusters;

    for (PfoVector::const_iterator iter = parentPfos.begin(), iterEnd = parentPfos.end(); iter != iterEnd; ++iter)
        LArPfoHelper::GetClusters(*iter, clusterHitType, parentClusters);

    for (PfoVector::const_iterator iter = daughterPfos.begin(), iterEnd = daughterPfos.end(); iter != iterEnd; ++iter)
        LArPfoHelper::GetClusters(*iter, clusterHitType, daughterClusters);

     // Select short parent clusters
    for (ClusterList::const_iterator cIter = parentClusters.begin(), cIterEnd = parentClusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;

        if (LArClusterHelper::GetLengthSquared(pCluster) > m_maxSeedClusterLength  * m_maxSeedClusterLength)
            continue;

        seedClusters.push_back(pCluster);
    }

    // Select all secondary clusters
    for (ClusterList::const_iterator cIter = daughterClusters.begin(), cIterEnd = daughterClusters.end(); cIter != cIterEnd; ++cIter)
    {
        seedClusters.push_back(*cIter);
    }

    // Select other possible delta rays
    for (ClusterVector::const_iterator cIter = inputClusters.begin(), cIterEnd = inputClusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;
 
        if (pCluster->GetNCaloHits() < m_minSeedClusterCaloHits)
            continue;

        const float parentDistance(parentClusters.empty() ? std::numeric_limits<float>::max() : 
            LArClusterHelper::GetClosestDistance(pCluster, parentClusters));
        const float daughterDistance(daughterClusters.empty() ? std::numeric_limits<float>::max() : 
            LArClusterHelper::GetClosestDistance(pCluster, daughterClusters));

        if (parentDistance < m_maxSeedClusterDisplacement && daughterDistance > m_maxSeedClusterDisplacement)
        {
            seedClusters.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayGrowingAlgorithm::GetPfos(const std::string inputPfoListName, PfoVector &pfoVector) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
        pfoVector.push_back(*pIter);

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{       
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentPfoListName", m_parentPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterPfoListName", m_daughterPfoListName));

    m_minCaloHitsPerCluster = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_minSeedClusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeedClusterCaloHits", m_minSeedClusterCaloHits));

    m_maxSeedClusterLength = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterLength", m_maxSeedClusterLength));

    m_maxSeedClusterDisplacement = 1.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedClusterDisplacement", m_maxSeedClusterDisplacement));

    return ClusterGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
