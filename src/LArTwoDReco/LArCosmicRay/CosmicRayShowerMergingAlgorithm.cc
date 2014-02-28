/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerMergingAlgorithm.cc
 *
 *  @brief  Implementation of the cosimic ray shower merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

void CosmicRayShowerMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMergingAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    // Separate clusters into seed and non-seed lists
    ClusterVector seedVector, nonSeedVector;
    this->SeparateClusters(clusterVector, seedVector, nonSeedVector);

    // Populate intermediate map of associations between clusters
    ClusterMergeMap intermediateMergeMap;
    this->FillClusterMergeMap(nonSeedVector, intermediateMergeMap);
    this->FillClusterMergeMap(seedVector, nonSeedVector, intermediateMergeMap);

    // Populate final map of associations between clusters
    this->FillClusterMergeMap(seedVector, seedVector, intermediateMergeMap, clusterMergeMap);
    this->FillClusterMergeMap(nonSeedVector, seedVector, intermediateMergeMap, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMergingAlgorithm::SeparateClusters(const ClusterVector &inputVector, ClusterVector &seedVector, ClusterVector &nonSeedVector) const
{
    ClusterVector cosmicVector, nonCosmicVector;

    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) > m_maxClusterLength * m_maxClusterLength)
            cosmicVector.push_back(pCluster);
        else if (pCluster->IsAvailable())
            nonCosmicVector.push_back(pCluster);
    }

    for (ClusterVector::const_iterator iterI = nonCosmicVector.begin(), iterEndI = nonCosmicVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pClusterI = *iterI;

        bool isSeed(false);

        for (ClusterVector::const_iterator iterJ = cosmicVector.begin(), iterEndJ = cosmicVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster* pClusterJ = *iterJ;

            if (LArClusterHelper::GetClosestDistance(pClusterI, pClusterJ) < m_maxSeedDisplacement)
            {
                isSeed = true;
                break;
            }
        }

        if (isSeed)
            seedVector.push_back(pClusterI);
        else
            nonSeedVector.push_back(pClusterI);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMergingAlgorithm::FillClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &outputMergeMap) const
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
                outputMergeMap[pClusterI].insert(pClusterJ);
                outputMergeMap[pClusterJ].insert(pClusterI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMergingAlgorithm::FillClusterMergeMap(const ClusterVector &firstVector, const ClusterVector &secondVector,
    ClusterMergeMap &outputMergeMap) const
{
    for (ClusterVector::const_iterator iterI = firstVector.begin(), iterEndI = firstVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = secondVector.begin(), iterEndJ = secondVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            if (this->IsAssociated(pClusterI, pClusterJ))
            {
                outputMergeMap[pClusterI].insert(pClusterJ);
                outputMergeMap[pClusterJ].insert(pClusterI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerMergingAlgorithm::FillClusterMergeMap(const ClusterVector &seedVector, const ClusterVector &vetoVector,
    const ClusterMergeMap &inputMergeMap, ClusterMergeMap &outputMergeMap) const
{
    for (ClusterVector::const_iterator iterI = seedVector.begin(), iterEndI = seedVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pClusterI = *iterI;

        ClusterList associatedList;

        this->CollectAssociatedClusters(pClusterI, inputMergeMap, associatedList);

        if (associatedList.empty())
            continue;

        bool isGoodMerge(true);

        for (ClusterVector::const_iterator iterJ = vetoVector.begin(), iterEndJ = vetoVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            if (associatedList.count(pClusterJ))
            {
                isGoodMerge = false;
                break;
            }
        }

        if (isGoodMerge)
        {

// ClusterList tempList;
// tempList.insert(pClusterI);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "SeedCluster", RED);
// PandoraMonitoringApi::VisualizeClusters(&associatedList, "AssociatedClusters", BLUE);
// PandoraMonitoringApi::ViewEvent();

            for (ClusterList::iterator iterJ = associatedList.begin(), iterEndJ = associatedList.end(); iterJ != iterEndJ; ++iterJ)
            {
                Cluster *pClusterJ = *iterJ;
                outputMergeMap[pClusterI].insert(pClusterJ);
                outputMergeMap[pClusterJ].insert(pClusterI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayShowerMergingAlgorithm::IsAssociated(const Cluster* const pClusterI, const Cluster* const pClusterJ) const
{
    if (LArClusterHelper::GetClosestDistance(pClusterI,pClusterJ) > m_maxNonSeedDisplacement)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLength = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    m_maxClusterLength = 15.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    m_maxSeedDisplacement = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeedDisplacement", m_maxSeedDisplacement));

    m_maxNonSeedDisplacement = 1.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxNonSeedDisplacement", m_maxNonSeedDisplacement));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
