/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the cosimic ray shower growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{


void CosmicRayShowerGrowingAlgorithm::GetListOfSeedClusters(const ClusterVector &inputClusters, ClusterVector &seedClusters) const
{
    if (inputClusters.empty())
        return;

    // Get hit type
    const Cluster* pFirstCluster = *(inputClusters.begin());
    const HitType clusterHitType(LArThreeDHelper::GetClusterHitType(pFirstCluster));

    // Select seed clusters for growing
    PfoVector primaryPfos, secondaryPfos;
    this->GetPfos(m_primaryPfoListName, primaryPfos);
    this->GetPfos(m_secondaryPfoListName, secondaryPfos);
    this->SelectPrimarySeeds(primaryPfos, clusterHitType, seedClusters);
    this->SelectSecondarySeeds(secondaryPfos, clusterHitType, seedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::SelectPrimarySeeds(const PfoVector &primaryPfos, const HitType clusterHitType,
    ClusterVector &seedClusters) const
{
    for (PfoVector::const_iterator pIter = primaryPfos.begin(), pIterEnd = primaryPfos.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;

        if (!pPfo->GetDaughterPfoList().empty())
            continue;

        ClusterList pfoClusterList;
        this->GetPfoClusters(pPfo, clusterHitType, pfoClusterList);

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pPfoCluster = *cIter;

            if (LArClusterHelper::GetLengthSquared(pPfoCluster) > m_maxPrimaryClusterLength  * m_maxPrimaryClusterLength)
                continue;

            seedClusters.push_back(pPfoCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayShowerGrowingAlgorithm::SelectSecondarySeeds(const PfoVector &secondaryPfos, const HitType clusterHitType,
    ClusterVector &seedClusters) const
{
    ClusterList pfoClusterList;
    this->GetPfoClusters(secondaryPfos, clusterHitType, pfoClusterList);

    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        seedClusters.push_back(*cIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxPrimaryClusterLength = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPrimaryClusterLength", m_maxPrimaryClusterLength));

    return CosmicRayGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
