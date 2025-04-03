/**
 *  @file   larpandoracontent/LArObjects/LArOverlapMatrix.cc
 *
 *  @brief  Implementation of the lar overlap matrix class.
 *
 *  $Log: $
 */

#include "Pandora/PandoraInputTypes.h"

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "Objects/Cluster.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArOverlapMatrix.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

template <typename T>
void OverlapMatrix<T>::GetUnambiguousElements(const bool ignoreUnavailable, ElementList &elementList) const
{
    for (typename TheMatrix::const_iterator iter1 = this->begin(), iter1End = this->end(); iter1 != iter1End; ++iter1)
    {
        ElementList tempElementList;
        ClusterList clusterList1, clusterList2;
        this->GetConnectedElements(iter1->first, ignoreUnavailable, tempElementList, clusterList1, clusterList2);

        const Cluster *pCluster1(nullptr), *pCluster2(nullptr);
        if (!this->DefaultAmbiguityFunction(clusterList1, clusterList2, pCluster1, pCluster2))
            continue;

        // ATTN With HIT_CUSTOM definitions, it is possible to navigate from different view 1 clusters to same combination
        if (iter1->first != pCluster1)
            continue;

        if (!pCluster1 || !pCluster2)
            continue;

        typename OverlapList::const_iterator iter2 = iter1->second.find(pCluster2);
        if (iter1->second.end() == iter2)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Element element(pCluster1, pCluster2, iter2->second);
        elementList.emplace_back(element);
    }

    std::sort(elementList.begin(), elementList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool OverlapMatrix<T>::DefaultAmbiguityFunction(
    const ClusterList &clusterList1, const ClusterList &clusterList2, const Cluster *&pCluster1, const Cluster *&pCluster2) const
{
    if ((1 != clusterList1.size()) || (1 != clusterList2.size()))
        return false;

    pCluster1 = *(clusterList1.begin());
    pCluster2 = *(clusterList2.begin());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::GetConnectedElements(
    const pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList, unsigned int &n1, unsigned int &n2) const
{
    ClusterList clusterList1, clusterList2;
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, clusterList1, clusterList2);
    n1 = clusterList1.size();
    n2 = clusterList2.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::GetSortedKeyClusters(ClusterVector &sortedKeyClusters) const
{
    for (typename TheMatrix::const_iterator iter1 = this->begin(), iter1End = this->end(); iter1 != iter1End; ++iter1)
        sortedKeyClusters.emplace_back(iter1->first);

    std::sort(sortedKeyClusters.begin(), sortedKeyClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::SetOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const OverlapResult &overlapResult)
{
    OverlapList &overlapList = m_overlapMatrix[pCluster1];
    typename OverlapList::const_iterator iter = overlapList.find(pCluster2);

    if (overlapList.end() != iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!overlapList.insert(typename OverlapList::value_type(pCluster2, overlapResult)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    ClusterList &navigation12(m_clusterNavigationMap12[pCluster1]);
    ClusterList &navigation21(m_clusterNavigationMap21[pCluster2]);

    if (navigation12.end() == std::find(navigation12.begin(), navigation12.end(), pCluster2))
        navigation12.emplace_back(pCluster2);
    if (navigation21.end() == std::find(navigation21.begin(), navigation21.end(), pCluster1))
        navigation21.emplace_back(pCluster1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::ReplaceOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const OverlapResult &overlapResult)
{
    typename TheMatrix::iterator iter1 = m_overlapMatrix.find(pCluster1);

    if (m_overlapMatrix.end() == iter1)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    typename OverlapList::iterator iter2 = iter1->second.find(pCluster2);

    if (iter1->second.end() == iter2)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    iter2->second = overlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::RemoveCluster(const pandora::Cluster *const pCluster)
{
    ClusterList additionalRemovals;

    if (m_clusterNavigationMap12.erase(pCluster) > 0)
    {
        typename TheMatrix::iterator iter = m_overlapMatrix.find(pCluster);

        if (m_overlapMatrix.end() != iter)
            m_overlapMatrix.erase(iter);

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMap21.begin(); navIter != m_clusterNavigationMap21.end();)
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = std::find(thisIter->second.begin(), thisIter->second.end(), pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.emplace_back(thisIter->first);
        }
    }

    if (m_clusterNavigationMap21.erase(pCluster) > 0)
    {
        for (typename TheMatrix::iterator iter1 = m_overlapMatrix.begin(), iter1End = m_overlapMatrix.end(); iter1 != iter1End; ++iter1)
        {
            typename OverlapList::iterator iter = iter1->second.find(pCluster);

            if (iter1->second.end() != iter)
                iter1->second.erase(iter);
        }

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMap12.begin(); navIter != m_clusterNavigationMap12.end();)
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = std::find(thisIter->second.begin(), thisIter->second.end(), pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.emplace_back(thisIter->first);
        }
    }

    additionalRemovals.sort(LArClusterHelper::SortByNHits);

    for (ClusterList::const_iterator iter = additionalRemovals.begin(), iterEnd = additionalRemovals.end(); iter != iterEnd; ++iter)
        this->RemoveCluster(*iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::GetConnectedElements(const Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
    ClusterList &clusterList1, ClusterList &clusterList2) const
{
    ClusterList localClusterList1, localClusterList2;
    this->ExploreConnections(pCluster, ignoreUnavailable, localClusterList1, localClusterList2);

    // ATTN Now need to check that all clusters received are from fully available matrix elements
    elementList.clear();
    clusterList1.clear();
    clusterList2.clear();

    for (typename TheMatrix::const_iterator iter1 = this->begin(), iter1End = this->end(); iter1 != iter1End; ++iter1)
    {
        if (localClusterList1.end() == std::find(localClusterList1.begin(), localClusterList1.end(), iter1->first))
            continue;

        for (typename OverlapList::const_iterator iter2 = iter1->second.begin(), iter2End = iter1->second.end(); iter2 != iter2End; ++iter2)
        {
            if (ignoreUnavailable && (!iter1->first->IsAvailable() || !iter2->first->IsAvailable()))
                continue;

            Element element(iter1->first, iter2->first, iter2->second);
            elementList.emplace_back(element);

            if (clusterList1.end() == std::find(clusterList1.begin(), clusterList1.end(), iter1->first))
                clusterList1.emplace_back(iter1->first);
            if (clusterList2.end() == std::find(clusterList2.begin(), clusterList2.end(), iter2->first))
                clusterList2.emplace_back(iter2->first);
        }
    }

    std::sort(elementList.begin(), elementList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapMatrix<T>::ExploreConnections(
    const Cluster *const pCluster, const bool ignoreUnavailable, ClusterList &clusterList1, ClusterList &clusterList2) const
{
    if (ignoreUnavailable && !pCluster->IsAvailable())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const bool clusterFromView1(
        !m_clusterNavigationMap12.empty() && (LArClusterHelper::GetClusterHitType(m_clusterNavigationMap12.begin()->first) == hitType));
    const bool clusterFromView2(
        !m_clusterNavigationMap21.empty() && (LArClusterHelper::GetClusterHitType(m_clusterNavigationMap21.begin()->first) == hitType));

    if (clusterFromView1 == clusterFromView2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList(clusterFromView1 ? clusterList1 : clusterList2);
    const ClusterNavigationMap &navigationMap(clusterFromView1 ? m_clusterNavigationMap12 : m_clusterNavigationMap21);

    if (clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pCluster))
        return;

    clusterList.emplace_back(pCluster);
    ClusterNavigationMap::const_iterator iter = navigationMap.find(pCluster);

    if (navigationMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterList::const_iterator cIter = iter->second.begin(), cIterEnd = iter->second.end(); cIter != cIterEnd; ++cIter)
        this->ExploreConnections(*cIter, ignoreUnavailable, clusterList1, clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class OverlapMatrix<float>;
template class OverlapMatrix<TwoViewDeltaRayOverlapResult>;
template class OverlapMatrix<TwoViewTransverseOverlapResult>;

} // namespace lar_content
