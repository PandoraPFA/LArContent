/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the n view matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArShowerOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingContainer.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingContainer.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
NViewMatchingAlgorithm<T>::NViewMatchingAlgorithm() :
    m_matchingContainer(this)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
NViewMatchingAlgorithm<T>::~NViewMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    m_matchingContainer.UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    m_matchingContainer.UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const std::string &NViewMatchingAlgorithm<T>::GetClusterListName(const HitType hitType) const
{
    return m_matchingContainer.GetClusterListName(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const pandora::ClusterList &NViewMatchingAlgorithm<T>::GetInputClusterList(const HitType hitType) const
{
    return m_matchingContainer.GetInputClusterList(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const pandora::ClusterList &NViewMatchingAlgorithm<T>::GetSelectedClusterList(const HitType hitType) const
{
    return m_matchingContainer.GetSelectedClusterList(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::SelectAllInputClusters()
{
    m_matchingContainer.SelectAllInputClusters();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::TidyUp()
{
    m_matchingContainer.TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::PerformMainLoop()
{
    m_matchingContainer.PerformMainLoop();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode NViewMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_matchingContainer.ReadSettings(xmlHandle));

    return MatchingBaseAlgorithm::ReadSettings(xmlHandle); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class NViewMatchingAlgorithm<TwoViewMatchingContainer<float> >;

template class NViewMatchingAlgorithm<ThreeViewMatchingContainer<float> >;
template class NViewMatchingAlgorithm<ThreeViewMatchingContainer<ShowerOverlapResult> >;

} // namespace lar_content
