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
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
NViewMatchingAlgorithm<T>::NViewMatchingAlgorithm() :
    m_matchingControl(this)
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
    m_matchingControl.UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    m_matchingControl.UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::string &NViewMatchingAlgorithm<T>::GetClusterListName(const HitType hitType) const
{
    return m_matchingControl.GetClusterListName(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const pandora::ClusterList &NViewMatchingAlgorithm<T>::GetInputClusterList(const HitType hitType) const
{
    return m_matchingControl.GetInputClusterList(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const pandora::ClusterList &NViewMatchingAlgorithm<T>::GetSelectedClusterList(const HitType hitType) const
{
    return m_matchingControl.GetSelectedClusterList(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::SelectAllInputClusters()
{
    m_matchingControl.SelectAllInputClusters();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::PrepareAllInputClusters()
{
    m_matchingControl.PrepareAllInputClusters();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::TidyUp()
{
    m_matchingControl.TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewMatchingAlgorithm<T>::PerformMainLoop()
{
    m_matchingControl.PerformMainLoop();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode NViewMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_matchingControl.ReadSettings(xmlHandle));

    return MatchingBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class NViewMatchingAlgorithm<TwoViewMatchingControl<float>>;
template class NViewMatchingAlgorithm<TwoViewMatchingControl<TwoViewDeltaRayOverlapResult>>;
template class NViewMatchingAlgorithm<TwoViewMatchingControl<TwoViewTransverseOverlapResult>>;
template class NViewMatchingAlgorithm<ThreeViewMatchingControl<float>>;
template class NViewMatchingAlgorithm<ThreeViewMatchingControl<ShowerOverlapResult>>;
template class NViewMatchingAlgorithm<ThreeViewMatchingControl<TransverseOverlapResult>>;
template class NViewMatchingAlgorithm<ThreeViewMatchingControl<LongitudinalOverlapResult>>;
template class NViewMatchingAlgorithm<ThreeViewMatchingControl<FragmentOverlapResult>>;
template class NViewMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult>>;

} // namespace lar_content
