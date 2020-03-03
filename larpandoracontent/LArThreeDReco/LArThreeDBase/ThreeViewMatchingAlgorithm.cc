/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the matching base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArOverlapTensor.h"
#include "larpandoracontent/LArObjects/LArShowerOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
ThreeViewMatchingAlgorithm<T>::ThreeViewMatchingAlgorithm() :
    m_pInputClusterListU(nullptr),
    m_pInputClusterListV(nullptr),
    m_pInputClusterListW(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
ThreeViewMatchingAlgorithm<T>::~ThreeViewMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingAlgorithm<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
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
                this->CalculateOverlapResult(pNewCluster, pCluster2, pCluster3);
            }
            else if (TPC_VIEW_V == hitType)
            {
                this->CalculateOverlapResult(pCluster2, pNewCluster, pCluster3);
            }
            else
            {
                this->CalculateOverlapResult(pCluster2, pCluster3, pNewCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
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

template<typename T>
const std::string &ThreeViewMatchingAlgorithm<T>::GetClusterListName(const HitType hitType) const
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

template<typename T>
const pandora::ClusterList &ThreeViewMatchingAlgorithm<T>::GetInputClusterList(const HitType hitType) const
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

template<typename T>
const pandora::ClusterList &ThreeViewMatchingAlgorithm<T>::GetSelectedClusterList(const HitType hitType) const
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
void ThreeViewMatchingAlgorithm<T>::SelectAllInputClusters()
{
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListNameU, m_pInputClusterListU));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListNameV, m_pInputClusterListV));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListNameW, m_pInputClusterListW));

    if (!m_pInputClusterListU || !m_pInputClusterListV || !m_pInputClusterListW)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ThreeViewMatchingAlgorithm: one or more input cluster lists unavailable." << std::endl;

        throw StatusCodeException(STATUS_CODE_SUCCESS);
    }

    this->SelectInputClusters(m_pInputClusterListU, m_clusterListU);
    this->SelectInputClusters(m_pInputClusterListV, m_clusterListV);
    this->SelectInputClusters(m_pInputClusterListW, m_clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewMatchingAlgorithm<T>::TidyUp()
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
void ThreeViewMatchingAlgorithm<T>::PerformMainLoop()
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
                this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode ThreeViewMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    return MatchingBaseAlgorithm::ReadSettings(xmlHandle);
}

template class ThreeViewMatchingAlgorithm<float>;
template class ThreeViewMatchingAlgorithm<TransverseOverlapResult>;
template class ThreeViewMatchingAlgorithm<LongitudinalOverlapResult>;
template class ThreeViewMatchingAlgorithm<FragmentOverlapResult>;
template class ThreeViewMatchingAlgorithm<ShowerOverlapResult>;

} // namespace lar_content
