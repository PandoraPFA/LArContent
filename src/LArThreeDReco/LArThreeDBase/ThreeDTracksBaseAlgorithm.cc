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
    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pNewCluster, m_slidingFitWindow, slidingFitResult);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pNewCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ThreeDBaseAlgorithm<T>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeDTracksBaseAlgorithm<T>::UpdateUponDeletion(Cluster *const pDeletedCluster)
{
    TwoDSlidingFitResultMap::iterator iter = m_slidingFitResultMap.find(pDeletedCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);

    ThreeDBaseAlgorithm<T>::UpdateUponDeletion(pDeletedCluster);
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
StatusCode ThreeDTracksBaseAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    return ThreeDBaseAlgorithm<T>::ReadSettings(xmlHandle);
}

template class ThreeDTracksBaseAlgorithm<float>;
template class ThreeDTracksBaseAlgorithm<TrackOverlapResult>;
template class ThreeDTracksBaseAlgorithm<TransverseOverlapResult>;
template class ThreeDTracksBaseAlgorithm<LongitudinalOverlapResult>;

} // namespace lar
