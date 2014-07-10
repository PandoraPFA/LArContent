/**
 *  @file   LArContent/src/LArObjects/LArOverlapTensor.cc
 * 
 *  @brief  Implementation of the lar overlap tensor class.
 * 
 *  $Log: $
 */

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraSettings.h"

#include "Pandora/PandoraInternal.h"
#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "Objects/Cluster.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"
#include "LArObjects/LArShowerOverlapResult.h"
#include "LArObjects/LArTrackOverlapResult.h"

using namespace pandora;

namespace lar
{

template <typename T>
void OverlapTensor<T>::GetUnambiguousElements(const bool ignoreUnavailable, AmbiguityFunction *pAmbiguityFunction, ElementList &elementList) const
{
    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
    {
        ElementList tempElementList;
        ClusterList clusterListU, clusterListV, clusterListW;
        this->GetConnectedElements(iterU->first, ignoreUnavailable, tempElementList, clusterListU, clusterListV, clusterListW);

        Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
        if (!(*pAmbiguityFunction)(clusterListU, clusterListV, clusterListW, pClusterU, pClusterV, pClusterW))
            continue;

        // ATTN: With custom definitions, it is possible to navigate from different U clusters to same combination
        if (iterU->first != pClusterU)
            continue;

        if ((NULL == pClusterU) || (NULL == pClusterV) || (NULL == pClusterW))
            continue;

        typename OverlapMatrix::const_iterator iterV = iterU->second.find(pClusterV);
        if (iterU->second.end() == iterV)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        typename OverlapList::const_iterator iterW = iterV->second.find(pClusterW);
        if (iterV->second.end() == iterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Element element(pClusterU, pClusterV, pClusterW, iterW->second);
        elementList.push_back(element);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool OverlapTensor<T>::DefaultAmbiguityFunction(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    if ((1 != clusterListU.size()) || (1 != clusterListV.size()) || (1 != clusterListW.size()))
        return false;

    pClusterU = *(clusterListU.begin());
    pClusterV = *(clusterListV.begin());
    pClusterW = *(clusterListW.begin());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
    unsigned int &nU, unsigned int &nV, unsigned int &nW) const
{
    ClusterList clusterListU, clusterListV, clusterListW;
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, clusterListU, clusterListV, clusterListW);
    nU = clusterListU.size(); nV = clusterListV.size(); nW = clusterListW.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
    pandora::Cluster *pClusterW, const OverlapResult &overlapResult)
{
    OverlapList &overlapList = m_overlapTensor[pClusterU][pClusterV];
    typename OverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() != iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!overlapList.insert(typename OverlapList::value_type(pClusterW, overlapResult)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    m_clusterNavigationMapUV[pClusterU].insert(pClusterV);
    m_clusterNavigationMapVW[pClusterV].insert(pClusterW);
    m_clusterNavigationMapWU[pClusterW].insert(pClusterU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::ReplaceOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
    pandora::Cluster *pClusterW, const OverlapResult &overlapResult)
{
    typename TheTensor::iterator iterU = m_overlapTensor.find(pClusterU);

    if (m_overlapTensor.end() == iterU)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    typename OverlapMatrix::iterator iterV = iterU->second.find(pClusterV);

    if (iterU->second.end() == iterV)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    typename OverlapList::iterator iterW = iterV->second.find(pClusterW);

    if (iterV->second.end() == iterW)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    iterW->second = overlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::RemoveCluster(pandora::Cluster *pCluster)
{
    ClusterList additionalRemovals;

    if (m_clusterNavigationMapUV.erase(pCluster) > 0)
    {
        typename TheTensor::iterator iter = m_overlapTensor.find(pCluster);

        if (m_overlapTensor.end() != iter)
            m_overlapTensor.erase(iter);

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMapWU.begin(); navIter != m_clusterNavigationMapWU.end(); )
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = thisIter->second.find(pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.insert(thisIter->first);
        }
    }

    if (m_clusterNavigationMapVW.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            typename OverlapMatrix::iterator iter = iterU->second.find(pCluster);

            if (iterU->second.end() != iter)
                iterU->second.erase(iter);
        }

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMapUV.begin(); navIter != m_clusterNavigationMapUV.end(); )
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = thisIter->second.find(pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.insert(thisIter->first);
        }
    }

    if (m_clusterNavigationMapWU.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            for (typename OverlapMatrix::iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
            {
                typename OverlapList::iterator iter = iterV->second.find(pCluster);

                if (iterV->second.end() != iter)
                    iterV->second.erase(iter);
            }
        }

        for (ClusterNavigationMap::iterator navIter = m_clusterNavigationMapVW.begin(); navIter != m_clusterNavigationMapVW.end(); )
        {
            ClusterNavigationMap::iterator thisIter = navIter++;
            ClusterList::iterator listIter = thisIter->second.find(pCluster);

            if (thisIter->second.end() != listIter)
                thisIter->second.erase(listIter);

            if (thisIter->second.empty())
                additionalRemovals.insert(thisIter->first);
        }
    }

    for (ClusterList::const_iterator iter = additionalRemovals.begin(), iterEnd = additionalRemovals.end(); iter != iterEnd; ++iter)
    {
        this->RemoveCluster(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedElements(Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
    ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    ClusterList localClusterListU, localClusterListV, localClusterListW;
    this->ExploreConnections(pCluster, ignoreUnavailable, localClusterListU, localClusterListV, localClusterListW);

    // ATTN Now need to check that all clusters received are from fully available tensor elements
    elementList.clear(); clusterListU.clear(); clusterListV.clear(); clusterListW.clear();

    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
    {
        if (0 == localClusterListU.count(iterU->first))
            continue;

        for (typename OverlapMatrix::const_iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
        {
            for (typename OverlapList::const_iterator iterW = iterV->second.begin(), iterWEnd = iterV->second.end(); iterW != iterWEnd; ++iterW)
            {
                if (ignoreUnavailable && (!iterU->first->IsAvailable() || !iterV->first->IsAvailable() || !iterW->first->IsAvailable()))
                    continue;

                Element element(iterU->first, iterV->first, iterW->first, iterW->second);
                elementList.push_back(element);

                clusterListU.insert(iterU->first);
                clusterListV.insert(iterV->first);
                clusterListW.insert(iterW->first);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::ExploreConnections(Cluster *const pCluster, const bool ignoreUnavailable, ClusterList &clusterListU,
    ClusterList &clusterListV, ClusterList &clusterListW) const
{
    if (ignoreUnavailable && !pCluster->IsAvailable())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
    const ClusterNavigationMap &navigationMap((TPC_VIEW_U == hitType) ? m_clusterNavigationMapUV : (TPC_VIEW_V == hitType) ? m_clusterNavigationMapVW : m_clusterNavigationMapWU);

    if (!clusterList.insert(pCluster).second)
        return;

    ClusterNavigationMap::const_iterator iter = navigationMap.find(pCluster);

    if (navigationMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterList::const_iterator cIter = iter->second.begin(), cIterEnd = iter->second.end(); cIter != cIterEnd; ++cIter)
        this->ExploreConnections(*cIter, ignoreUnavailable, clusterListU, clusterListV, clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class OverlapTensor<float>;
template class OverlapTensor<TransverseOverlapResult>;
template class OverlapTensor<LongitudinalOverlapResult>;
template class OverlapTensor<FragmentOverlapResult>;
template class OverlapTensor<ShowerOverlapResult>;

} // namespace lar
