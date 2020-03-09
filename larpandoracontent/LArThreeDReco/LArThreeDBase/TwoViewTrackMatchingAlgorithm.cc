/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the two view track matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewTrackMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template<typename T>
TwoViewTrackMatchingAlgorithm<T>::TwoViewTrackMatchingAlgorithm() :
    m_slidingFitWindow(20),
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TwoViewTrackMatchingAlgorithm<T>::~TwoViewTrackMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const TwoDSlidingFitResult &TwoViewTrackMatchingAlgorithm<T>::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    try
    {
        this->AddToSlidingFitCache(pNewCluster);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;

        return;
    }

    TwoViewMatchingAlgorithm<T>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    TwoViewMatchingAlgorithm<T>::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewTrackMatchingAlgorithm<T>::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (const Cluster *const pCluster : *pInputClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::SetPfoParticleId(PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    pfoParameters.m_particleId = MU_MINUS; // Track
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::RemoveFromSlidingFitCache(const Cluster *const pCluster)
{
    TwoDSlidingFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::PreparationStep(ClusterList &clusterList)
{
    for (ClusterList::iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; )
    {
        try
        {
            this->AddToSlidingFitCache(*iter);
            ++iter;
        }
        catch (const StatusCodeException &statusCodeException)
        {
            clusterList.erase(iter++);

            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::PreparationStep()
{
    this->PreparationStep(this->m_clusterList1);
    this->PreparationStep(this->m_clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void TwoViewTrackMatchingAlgorithm<T>::TidyUp()
{
    m_slidingFitResultMap.clear();
    return TwoViewMatchingAlgorithm<T>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode TwoViewTrackMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    return TwoViewMatchingAlgorithm<T>::ReadSettings(xmlHandle);
}

template class TwoViewTrackMatchingAlgorithm<float>;

} // namespace lar_content
