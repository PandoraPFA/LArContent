/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/ClusterMopUpBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster mop up algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/ClusterMopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterMopUpBaseAlgorithm::ClusterMopUpBaseAlgorithm() :
    m_excludePfosContainingTracks(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMopUpBaseAlgorithm::Run()
{
    ClusterList pfoClusterListU, pfoClusterListV, pfoClusterListW;
    this->GetPfoClusterLists(pfoClusterListU, pfoClusterListV, pfoClusterListW);

    ClusterList remnantClusterListU, remnantClusterListV, remnantClusterListW;
    this->GetRemnantClusterLists(remnantClusterListU, remnantClusterListV, remnantClusterListW);

    this->ClusterMopUp(pfoClusterListU, remnantClusterListU);
    this->ClusterMopUp(pfoClusterListV, remnantClusterListV);
    this->ClusterMopUp(pfoClusterListW, remnantClusterListW);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpBaseAlgorithm::GetPfoClusterLists(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (StringVector::const_iterator sIter = m_pfoListNames.begin(), sIterEnd = m_pfoListNames.end(); sIter != sIterEnd; ++sIter)
    {
        const PfoList *pPfoList = NULL;
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *sIter, pPfoList))
            continue;

        for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
        {
            const ParticleFlowObject *const pPfo = *pIter;

            if (m_excludePfosContainingTracks && LArPfoHelper::IsTrack(pPfo))
                continue;

            this->GetClusterLists(pPfo->GetClusterList(), false, clusterListU, clusterListV, clusterListW);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpBaseAlgorithm::GetRemnantClusterLists(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (StringVector::const_iterator sIter = m_remnantClusterListNames.begin(), sIterEnd = m_remnantClusterListNames.end(); sIter != sIterEnd; ++sIter)
    {
        const ClusterList *pClusterList(NULL);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *sIter, pClusterList))
            continue;

        this->GetClusterLists(*pClusterList, true, clusterListU, clusterListV, clusterListW);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpBaseAlgorithm::GetClusterLists(const ClusterList &inputClusterList, const bool availabilityFlag, ClusterList &clusterListU,
    ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (ClusterList::const_iterator cIter = inputClusterList.begin(), cIterEnd = inputClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster(*cIter);

        if (availabilityFlag != pCluster->IsAvailable())
            continue;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            continue;

        ClusterList &target((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
        target.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpBaseAlgorithm::MakeClusterMerges(const ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterVector sortedRemnantClusters;
    for (const auto &remnantMapEntry : clusterAssociationMap) sortedRemnantClusters.push_back(remnantMapEntry.first);
    std::sort(sortedRemnantClusters.begin(), sortedRemnantClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pRemnantCluster : sortedRemnantClusters)
    {
        const AssociationDetails &associationDetails(clusterAssociationMap.at(pRemnantCluster));
        const Cluster *pBestPfoCluster(nullptr);
        float bestFigureOfMerit(-std::numeric_limits<float>::max());

        ClusterVector sortedPfoClusters;
        for (const auto &pfoMapEntry : associationDetails) sortedPfoClusters.push_back(pfoMapEntry.first);
        std::sort(sortedPfoClusters.begin(), sortedPfoClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pPfoCluster : sortedPfoClusters)
        {
            const float figureOfMerit(associationDetails.at(pPfoCluster));

            if (figureOfMerit > bestFigureOfMerit)
            {
                pBestPfoCluster = pPfoCluster;
                bestFigureOfMerit = figureOfMerit;
            }
        }

        if (!pBestPfoCluster)
            continue;

        const std::string listNameP(this->GetListName(pBestPfoCluster));
        const std::string listNameR(this->GetListName(pRemnantCluster));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pBestPfoCluster, pRemnantCluster,
            listNameP, listNameR));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::string ClusterMopUpBaseAlgorithm::GetListName(const T *const pT) const
{
    std::string currentListName;
    const std::MANAGED_CONTAINER<const T*> *pCurrentList(nullptr);
    (void) PandoraContentApi::GetCurrentList(*this, pCurrentList, currentListName);

    if (pCurrentList && (pCurrentList->count(pT)))
        return currentListName;

    for (const std::string &listName : m_remnantClusterListNames)
    {
        const std::MANAGED_CONTAINER<const T*> *pList(nullptr);
        (void) PandoraContentApi::GetList(*this, listName, pList);

        if (pList && (pList->count(pT)))
            return listName;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMopUpBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "RemnantClusterListNames", m_remnantClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExcludePfosContainingTracks", m_excludePfosContainingTracks));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
