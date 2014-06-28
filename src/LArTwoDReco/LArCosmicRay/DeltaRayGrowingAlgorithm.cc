/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{


void DeltaRayGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    if (inputClusters.empty())
        return;

    // Get hit type
    const Cluster* pFirstCluster = *(inputClusters.begin());
    const HitType clusterHitType(LArThreeDHelper::GetClusterHitType(pFirstCluster));

    // Select seed clusters for growing
    PfoVector primaryPfos;
    this->GetPfos(m_primaryPfoListName, primaryPfos);
    this->SelectClusterSeeds(inputClusters, primaryPfos, clusterHitType, seedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayGrowingAlgorithm::SelectClusterSeeds(const ClusterVector &inputClusters, const PfoVector &primaryPfos,
    const HitType clusterHitType, ClusterVector &seedClusters) const
{
    ClusterList pfoClusterList;
    this->GetPfoClusters(primaryPfos, clusterHitType, pfoClusterList);

    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (clusterHitType != LArThreeDHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        bool isSeed(false);

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *pPfoCluster = *cIter;

            if (LArClusterHelper::GetClosestDistance(pCluster, pPfoCluster) < m_maxPrimaryClusterDisplacement)
            {
                isSeed = true;
                break;
            }
        }

        if (isSeed)
            seedClusters.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{    
    m_maxPrimaryClusterDisplacement = 1.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPrimaryClusterDisplacement", m_maxPrimaryClusterDisplacement));

    return CosmicRayGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
