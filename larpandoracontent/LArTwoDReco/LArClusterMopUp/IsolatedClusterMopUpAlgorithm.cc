/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/IsolatedClusterMopUpAlgorithm.cc
 *
 *  @brief  Implementation of the isolated cluster mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/IsolatedClusterMopUpAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

IsolatedClusterMopUpAlgorithm::IsolatedClusterMopUpAlgorithm() :
    m_maxCaloHitsInCluster(20), m_maxHitClusterDistance(5.f), m_addHitsAsIsolated(true)
{
    // ATTN Default value differs from base class
    m_excludePfosContainingTracks = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedClusterMopUpAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters) const
{
    CaloHitList caloHitList;
    this->DissolveClustersToHits(remnantClusters, caloHitList);

    // ATTN remnantClusters now contains dangling pointers
    CaloHitToClusterMap caloHitToClusterMap;
    this->GetCaloHitToClusterMap(caloHitList, pfoClusters, caloHitToClusterMap);

    CaloHitList sortedCaloHitList;
    for (const auto &mapEntry : caloHitToClusterMap)
        sortedCaloHitList.push_back(mapEntry.first);
    sortedCaloHitList.sort(LArClusterHelper::SortHitsByPosition);

    for (const CaloHit *pCaloHit : sortedCaloHitList)
    {
        if (m_addHitsAsIsolated)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, caloHitToClusterMap.at(pCaloHit), pCaloHit));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, caloHitToClusterMap.at(pCaloHit), pCaloHit));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedClusterMopUpAlgorithm::DissolveClustersToHits(const ClusterList &clusterList, CaloHitList &caloHitList) const
{
    for (const Cluster *const pRemnantCluster : clusterList)
    {
        if (pRemnantCluster->GetNCaloHits() < m_maxCaloHitsInCluster)
        {
            const std::string listNameR(this->GetListName(pRemnantCluster));
            pRemnantCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pRemnantCluster, listNameR));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedClusterMopUpAlgorithm::GetCaloHitToClusterMap(
    const CaloHitList &caloHitList, const ClusterList &clusterList, CaloHitToClusterMap &caloHitToClusterMap) const
{
    CaloHitList allCaloHits;
    CaloHitToClusterMap hitToParentClusterMap;

    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);
        allCaloHits.insert(allCaloHits.end(), daughterHits.begin(), daughterHits.end());

        for (const CaloHit *const pCaloHit : daughterHits)
            (void)hitToParentClusterMap.insert(CaloHitToClusterMap::value_type(pCaloHit, pCluster));
    }

    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    HitType view{TPC_3D};
    if (!caloHitList.empty())
        view = caloHitList.front()->GetHitType();
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
    const float maxHitClusterDistanceAdjusted{ratio * m_maxHitClusterDistance};

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const HitKDNode2D *pResultHit(nullptr);
        float resultDistance(std::numeric_limits<float>::max());
        const HitKDNode2D targetHit(pCaloHit, pCaloHit->GetPositionVector().GetX(), pCaloHit->GetPositionVector().GetZ());
        kdTree.findNearestNeighbour(targetHit, pResultHit, resultDistance);

        if (pResultHit && (resultDistance < maxHitClusterDistanceAdjusted))
            (void)caloHitToClusterMap.insert(CaloHitToClusterMap::value_type(pCaloHit, hitToParentClusterMap.at(pResultHit->data)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCaloHitsInCluster", m_maxCaloHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxHitClusterDistance", m_maxHitClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AddHitsAsIsolated", m_addHitsAsIsolated));

    return ClusterMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
