/**
 *  @file   LArContent/src/LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional tracks tracks algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

template<typename T>
const TwoDSlidingFitResult &ThreeDTracksBaseAlgorithm<T>::GetCachedSlidingFitResult(Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::UpdateForNewCluster(Cluster *const pNewCluster)
{
    this->AddToSlidingFitCache(pNewCluster);
    ThreeDBaseAlgorithm<T>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::UpdateUponDeletion(Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    ThreeDBaseAlgorithm<T>::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDTracksBaseAlgorithm<T>::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::PreparationStep()
{
    ClusterList allClustersList;
    allClustersList.insert(this->m_clusterListU.begin(), this->m_clusterListU.end());
    allClustersList.insert(this->m_clusterListV.begin(), this->m_clusterListV.end());
    allClustersList.insert(this->m_clusterListW.begin(), this->m_clusterListW.end());

    for (ClusterList::const_iterator iter = allClustersList.begin(), iterEnd = allClustersList.end(); iter != iterEnd; ++iter)
    {
        TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(*iter, m_slidingFitWindow, slidingFitResult);

        if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::TidyUp()
{
    m_slidingFitResultMap.clear();
    return ThreeDBaseAlgorithm<T>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::AddToSlidingFitCache(Cluster *const pCluster)
{
    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, m_slidingFitWindow, slidingFitResult);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::RemoveFromSlidingFitCache(Cluster *const pCluster)
{
    TwoDSlidingFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode ThreeDTracksBaseAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_minClusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    return ThreeDBaseAlgorithm<T>::ReadSettings(xmlHandle);
}

template class ThreeDTracksBaseAlgorithm<float>;
template class ThreeDTracksBaseAlgorithm<TransverseOverlapResult>;
template class ThreeDTracksBaseAlgorithm<LongitudinalOverlapResult>;
template class ThreeDTracksBaseAlgorithm<FragmentOverlapResult>;

} // namespace lar
