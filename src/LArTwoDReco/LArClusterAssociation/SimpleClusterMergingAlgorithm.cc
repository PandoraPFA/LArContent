/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the simple cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SimpleClusterMergingAlgorithm::SimpleClusterMergingAlgorithm() :
    m_minCaloHitsPerCluster(10),
    m_maxClusterSeparation(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterMergingAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            if (this->IsAssociated(pClusterI, pClusterJ))
            {
                clusterMergeMap[pClusterI].insert(pClusterJ);
                clusterMergeMap[pClusterJ].insert(pClusterI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SimpleClusterMergingAlgorithm::IsAssociated(const Cluster* const pClusterI, const Cluster* const pClusterJ) const
{
    if (LArClusterHelper::GetClosestDistance(pClusterI,pClusterJ) > m_maxClusterSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterSeparation", m_maxClusterSeparation));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
