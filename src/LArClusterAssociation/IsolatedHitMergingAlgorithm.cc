/**
 *  @file   LArContent/src/LArClusterAssociation/IsolatedHitMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the isolated hit merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/IsolatedHitMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode IsolatedHitMergingAlgorithm::Run()
{
    // TODO - remove this test code. Later maybe remove all non-seed clusters, having promoted important remnants
    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_INITIALIZED != statusCode))
        return statusCode;

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    ClusterList clustersToDelete;
    for (ClusterList::const_iterator iter = pNonSeedClusterList->begin(), iterEnd = pNonSeedClusterList->end(); iter != iterEnd; ++iter)    
    {
         if ((*iter)->GetNCaloHits() < 10)
            clustersToDelete.insert(*iter);
    }

    // Delete the current list of non-seed clusters to free the underlying hits
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeleteClusters(*this, clustersToDelete, m_nonSeedClusterListName));

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), Cluster::SortByInnerLayer);

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentCaloHitList(*this, pCaloHitList));

    CaloHitVector caloHitVector(pCaloHitList->begin(), pCaloHitList->end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), IsolatedHitMergingAlgorithm::SortByLayer);

    for (CaloHitList::const_iterator iterI = pCaloHitList->begin(), iterIEnd = pCaloHitList->end(); iterI != iterIEnd; ++iterI)
    {
        CaloHit *pCaloHit = *iterI;

        if (!PandoraContentApi::IsCaloHitAvailable(*this, pCaloHit))
            continue;

        Cluster *pBestHostCluster(NULL);
        float bestHostClusterEnergy(0.);
        float minDistance(10.f);

        // Find most appropriate cluster for this isolated hit
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
// ClusterList pTempList;
// CaloHitList dTempList;
// pTempList.insert(pBestHostCluster);
// dTempList.insert(pCaloHit);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&pTempList, "ISOParent", RED);
// PandoraMonitoringApi::VisualizeCaloHits(&dTempList, "ISODaughter", BLUE);
// PandoraMonitoringApi::ViewEvent();
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedCaloHitToCluster(*this, pBestHostCluster, pCaloHit));
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
