/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.cc
 *
 *  @brief  Implementation of the three view matching control class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArShowerOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
ThreeViewMatchingControl<T>::ThreeViewMatchingControl(MatchingBaseAlgorithm *const pAlgorithm) :
    NViewMatchingControl(pAlgorithm), m_pInputClusterListU(nullptr), m_pInputClusterListV(nullptr), m_pInputClusterListW(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
ThreeViewMatchingControl<T>::~ThreeViewMatchingControl()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
typename ThreeViewMatchingControl<T>::TensorType &ThreeViewMatchingControl<T>::GetOverlapTensor()
{
    return m_overlapTensor;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingControl<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pNewCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((TPC_VIEW_U == hitType) ? m_clusterListU : (TPC_VIEW_V == hitType) ? m_clusterListV : m_clusterListW);

    if (clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pNewCluster))
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    clusterList.push_back(pNewCluster);

    const ClusterList &clusterList2((TPC_VIEW_U == hitType) ? m_clusterListV : m_clusterListU);
    const ClusterList &clusterList3((TPC_VIEW_W == hitType) ? m_clusterListV : m_clusterListW);

    ClusterVector clusterVector2(clusterList2.begin(), clusterList2.end());
    ClusterVector clusterVector3(clusterList3.begin(), clusterList3.end());
    std::sort(clusterVector2.begin(), clusterVector2.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVector3.begin(), clusterVector3.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster2 : clusterVector2)
    {
        for (const Cluster *const pCluster3 : clusterVector3)
        {
            if (TPC_VIEW_U == hitType)
            {
                m_pAlgorithm->CalculateOverlapResult(pNewCluster, pCluster2, pCluster3);
            }
            else if (TPC_VIEW_V == hitType)
            {
                m_pAlgorithm->CalculateOverlapResult(pCluster2, pNewCluster, pCluster3);
            }
            else
            {
                m_pAlgorithm->CalculateOverlapResult(pCluster2, pCluster3, pNewCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingControl<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    ClusterList::iterator iterU = std::find(m_clusterListU.begin(), m_clusterListU.end(), pDeletedCluster);
    ClusterList::iterator iterV = std::find(m_clusterListV.begin(), m_clusterListV.end(), pDeletedCluster);
    ClusterList::iterator iterW = std::find(m_clusterListW.begin(), m_clusterListW.end(), pDeletedCluster);

    if (m_clusterListU.end() != iterU)
        m_clusterListU.erase(iterU);

    if (m_clusterListV.end() != iterV)
        m_clusterListV.erase(iterV);

    if (m_clusterListW.end() != iterW)
        m_clusterListW.erase(iterW);

    m_overlapTensor.RemoveCluster(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::string &ThreeViewMatchingControl<T>::GetClusterListName(const HitType hitType) const
{
    if (TPC_VIEW_U == hitType)
        return m_inputClusterListNameU;

    if (TPC_VIEW_V == hitType)
        return m_inputClusterListNameV;

    if (TPC_VIEW_W == hitType)
        return m_inputClusterListNameW;

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const pandora::ClusterList &ThreeViewMatchingControl<T>::GetInputClusterList(const HitType hitType) const
{
    if ((TPC_VIEW_U == hitType) && m_pInputClusterListU)
        return (*m_pInputClusterListU);

    if ((TPC_VIEW_V == hitType) && m_pInputClusterListV)
        return (*m_pInputClusterListV);

    if ((TPC_VIEW_W == hitType) && m_pInputClusterListW)
        return (*m_pInputClusterListW);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const pandora::ClusterList &ThreeViewMatchingControl<T>::GetSelectedClusterList(const HitType hitType) const
{
    if (TPC_VIEW_U == hitType)
        return m_clusterListU;

    if (TPC_VIEW_V == hitType)
        return m_clusterListV;

    if (TPC_VIEW_W == hitType)
        return m_clusterListW;

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingControl<T>::SelectAllInputClusters()
{
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
        PandoraContentApi::GetList(*m_pAlgorithm, m_inputClusterListNameU, m_pInputClusterListU));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
        PandoraContentApi::GetList(*m_pAlgorithm, m_inputClusterListNameV, m_pInputClusterListV));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
        PandoraContentApi::GetList(*m_pAlgorithm, m_inputClusterListNameW, m_pInputClusterListW));

    if (!m_pInputClusterListU || !m_pInputClusterListV || !m_pInputClusterListW)
    {
        if (PandoraContentApi::GetSettings(*m_pAlgorithm)->ShouldDisplayAlgorithmInfo())
            std::cout << "ThreeViewMatchingControl: one or more input cluster lists unavailable." << std::endl;

        throw StatusCodeException(STATUS_CODE_SUCCESS);
    }

    m_pAlgorithm->SelectInputClusters(m_pInputClusterListU, m_clusterListU);
    m_pAlgorithm->SelectInputClusters(m_pInputClusterListV, m_clusterListV);
    m_pAlgorithm->SelectInputClusters(m_pInputClusterListW, m_clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingControl<T>::PrepareAllInputClusters()
{
    m_pAlgorithm->PrepareInputClusters(m_clusterListU);
    m_pAlgorithm->PrepareInputClusters(m_clusterListV);
    m_pAlgorithm->PrepareInputClusters(m_clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingControl<T>::TidyUp()
{
    m_overlapTensor.Clear();

    m_pInputClusterListU = nullptr;
    m_pInputClusterListV = nullptr;
    m_pInputClusterListW = nullptr;

    m_clusterListU.clear();
    m_clusterListV.clear();
    m_clusterListW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingControl<T>::PerformMainLoop()
{
    ClusterVector clusterVectorU(m_clusterListU.begin(), m_clusterListU.end());
    ClusterVector clusterVectorV(m_clusterListV.begin(), m_clusterListV.end());
    ClusterVector clusterVectorW(m_clusterListW.begin(), m_clusterListW.end());
    std::sort(clusterVectorU.begin(), clusterVectorU.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVectorV.begin(), clusterVectorV.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVectorW.begin(), clusterVectorW.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pClusterU : clusterVectorU)
    {
        for (const Cluster *const pClusterV : clusterVectorV)
        {
            for (const Cluster *const pClusterW : clusterVectorW)
                m_pAlgorithm->CalculateOverlapResult(pClusterU, pClusterV, pClusterW);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode ThreeViewMatchingControl<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    return STATUS_CODE_SUCCESS;
}

template class ThreeViewMatchingControl<float>;
template class ThreeViewMatchingControl<TransverseOverlapResult>;
template class ThreeViewMatchingControl<LongitudinalOverlapResult>;
template class ThreeViewMatchingControl<FragmentOverlapResult>;
template class ThreeViewMatchingControl<ShowerOverlapResult>;
template class ThreeViewMatchingControl<DeltaRayOverlapResult>;

} // namespace lar_content
