/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the isolated hit merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode IsolatedHitMergingAlgorithm::Run()
{
    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_seedClusterListName, pSeedClusterList));

    // Delete small clusters from the current list of non-seeds to free the underlying hits
    ClusterList clustersToDelete;
    CaloHitList caloHitList;

    for (ClusterList::const_iterator iter = pNonSeedClusterList->begin(), iterEnd = pNonSeedClusterList->end(); iter != iterEnd; ++iter)    
    {
        Cluster *pCluster = *iter;

        if (pCluster->GetNCaloHits() < 10)
        {
            clustersToDelete.insert(pCluster);
            pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clustersToDelete, m_nonSeedClusterListName));

    // Reassign available hits to be most appropriate seed clusters
    ClusterVector clusterVector(pSeedClusterList->begin(), pSeedClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), IsolatedHitMergingAlgorithm::SortByLayer);

    for (CaloHitVector::const_iterator iterI = caloHitVector.begin(), iterIEnd = caloHitVector.end(); iterI != iterIEnd; ++iterI)
    {
        CaloHit *pCaloHit = *iterI;

        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        Cluster *pBestHostCluster(NULL);
        float bestHostClusterEnergy(0.), minDistance(m_maxHitClusterDistance);

        for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            Cluster *pCluster = *iterJ;

            if (NULL == pCluster)
                continue;

            const float distance(this->GetDistanceToHit(pCluster, pCaloHit));
            const float hostClusterEnergy(pCluster->GetHadronicEnergy());

            if ((distance < minDistance) || ((distance == minDistance) && (hostClusterEnergy > bestHostClusterEnergy)))
            {
                minDistance = distance;
                pBestHostCluster = pCluster;
                bestHostClusterEnergy = hostClusterEnergy;
            }
        }

        if (NULL != pBestHostCluster)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, pBestHostCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float IsolatedHitMergingAlgorithm::GetDistanceToHit(const Cluster *const pCluster, const CaloHit *const pCaloHit) const
{
    // Apply simple preselection using cosine of opening angle between the hit and cluster directions
    if (pCaloHit->GetExpectedDirection().GetCosOpeningAngle(pCluster->GetInitialDirection()) < 0.f)
        return std::numeric_limits<float>::max();

    float minDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        const float distanceSquared((pCluster->GetCentroid(iter->first) - hitPosition).GetMagnitudeSquared());

        if (distanceSquared < minDistanceSquared)
            minDistanceSquared = distanceSquared;
    }

    if (minDistanceSquared < std::numeric_limits<float>::max())
        return std::sqrt(minDistanceSquared);

    return std::numeric_limits<float>::max();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool IsolatedHitMergingAlgorithm::SortByLayer(const CaloHit *const pLhs, const CaloHit *const pRhs)
{
    const unsigned int layerLhs(pLhs->GetPseudoLayer());
    const unsigned int layerRhs(pRhs->GetPseudoLayer());

    if (layerLhs != layerRhs)
        return (layerLhs < layerRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    m_maxHitClusterDistance = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxHitClusterDistance", m_maxHitClusterDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
