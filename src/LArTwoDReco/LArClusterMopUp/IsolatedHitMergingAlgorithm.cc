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

namespace lar_content
{

void IsolatedHitMergingAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters,
    const ClusterToListNameMap &clusterToListNameMap) const
{
    CaloHitList caloHitList;
    this->DissolveClustersToHits(remnantClusters, clusterToListNameMap, caloHitList);

    CaloHitToClusterMap caloHitToClusterMap;
    this->GetCaloHitToClusterMap(caloHitList, pfoClusters, caloHitToClusterMap);

    for (CaloHitToClusterMap::const_iterator iter = caloHitToClusterMap.begin(), iterEnd = caloHitToClusterMap.end(); iter != iterEnd; ++iter)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, iter->second, iter->first));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedHitMergingAlgorithm::DissolveClustersToHits(const ClusterList &clusterList, const ClusterToListNameMap &clusterToListNameMap,
    CaloHitList &caloHitList) const
{
    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pRemnantCluster(*iter);

        if (pRemnantCluster->GetNCaloHits() < m_maxCaloHitsInCluster)
        {
            ClusterToListNameMap::const_iterator nameIter = clusterToListNameMap.find(pRemnantCluster);

            if (clusterToListNameMap.end() == nameIter)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            pRemnantCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pRemnantCluster, nameIter->second));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedHitMergingAlgorithm::GetCaloHitToClusterMap(const CaloHitList &caloHitList, const ClusterList &clusterList, CaloHitToClusterMap &caloHitToClusterMap) const
{
    for (CaloHitList::const_iterator hIter = caloHitList.begin(), hIterEnd = caloHitList.end(); hIter != hIterEnd; ++hIter)
    {
        CaloHit *pCaloHit(*hIter);

        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        float bestDistance(m_maxHitClusterDistance), tieBreakerBestEnergy(0.f);

        for (ClusterList::const_iterator pIter = clusterList.begin(), pIterEnd = clusterList.end(); pIter != pIterEnd; ++pIter)
        {
            Cluster *pPfoCluster(*pIter);
            const float clusterDistance(this->GetDistanceToHit(pPfoCluster, pCaloHit));
            const float clusterEnergy(pPfoCluster->GetElectromagneticEnergy());

            if ((clusterDistance < bestDistance) || ((clusterDistance == bestDistance) && (clusterEnergy > tieBreakerBestEnergy)))
            {
                bestDistance = clusterDistance;
                tieBreakerBestEnergy = clusterEnergy;
                caloHitToClusterMap[pCaloHit] = pPfoCluster;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float IsolatedHitMergingAlgorithm::GetDistanceToHit(const Cluster *const pCluster, const CaloHit *const pCaloHit) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const int hitLayer(static_cast<int>(pCaloHit->GetPseudoLayer()));
    const int clusterMinLayer(static_cast<int>(orderedCaloHitList.begin()->first));
    const int clusterMaxLayer(static_cast<int>(orderedCaloHitList.rbegin()->first));

    const int startLayer(std::max(clusterMinLayer, hitLayer - m_hitLayerSearchWindow));
    const int endLayer(std::min(clusterMaxLayer, hitLayer + m_hitLayerSearchWindow));

    float minDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

    for (int iLayer = startLayer; iLayer <= endLayer; ++iLayer)
    {
        OrderedCaloHitList::const_iterator iter = orderedCaloHitList.find(iLayer);

        if (orderedCaloHitList.end() == iter)
            continue;

        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const float distanceSquared(((*hIter)->GetPositionVector() - hitPosition).GetMagnitudeSquared());

            if (distanceSquared < minDistanceSquared)
                minDistanceSquared = distanceSquared;
        }
    }

    if (minDistanceSquared < std::numeric_limits<float>::max())
        return std::sqrt(minDistanceSquared);

    return std::numeric_limits<float>::max();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxCaloHitsInCluster = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxCaloHitsInCluster", m_maxCaloHitsInCluster));

    m_hitLayerSearchWindow = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "HitLayerSearchWindow", m_hitLayerSearchWindow));

    m_maxHitClusterDistance = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxHitClusterDistance", m_maxHitClusterDistance));

    return ClusterMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
