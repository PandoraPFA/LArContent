/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the n view track matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
NViewTrackMatchingAlgorithm<T>::NViewTrackMatchingAlgorithm() :
    m_slidingFitWindow(20), m_minClusterCaloHits(5), m_minClusterLengthSquared(3.f * 3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
NViewTrackMatchingAlgorithm<T>::~NViewTrackMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const TwoDSlidingFitResult &NViewTrackMatchingAlgorithm<T>::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool NViewTrackMatchingAlgorithm<T>::MakeClusterSplits(const SplitPositionMap &splitPositionMap)
{
    bool changesMade(false);

    ClusterList splitClusters;
    for (const auto &mapEntry : splitPositionMap)
        splitClusters.push_back(mapEntry.first);
    splitClusters.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *pCurrentCluster : splitClusters)
    {
        CartesianPointVector splitPositions(splitPositionMap.at(pCurrentCluster));
        std::sort(splitPositions.begin(), splitPositions.end(), NViewTrackMatchingAlgorithm<T>::SortSplitPositions);

        const HitType hitType(LArClusterHelper::GetClusterHitType(pCurrentCluster));
        const std::string &clusterListName(this->GetClusterListName(hitType));

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        for (const CartesianVector &splitPosition : splitPositions)
        {
            const Cluster *pLowXCluster(nullptr), *pHighXCluster(nullptr);
            this->UpdateUponDeletion(pCurrentCluster);

            if (this->MakeClusterSplit(splitPosition, pCurrentCluster, pLowXCluster, pHighXCluster))
            {
                changesMade = true;
                this->UpdateForNewCluster(pLowXCluster);
                this->UpdateForNewCluster(pHighXCluster);
                pCurrentCluster = pHighXCluster;
            }
            else
            {
                this->UpdateForNewCluster(pCurrentCluster);
            }
        }
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool NViewTrackMatchingAlgorithm<T>::MakeClusterSplit(
    const CartesianVector &splitPosition, const Cluster *&pCurrentCluster, const Cluster *&pLowXCluster, const Cluster *&pHighXCluster) const
{
    CartesianVector lowXEnd(0.f, 0.f, 0.f), highXEnd(0.f, 0.f, 0.f);

    try
    {
        LArPointingCluster pointingCluster(pCurrentCluster);
        const bool innerIsLowX(pointingCluster.GetInnerVertex().GetPosition().GetX() < pointingCluster.GetOuterVertex().GetPosition().GetX());
        lowXEnd = (innerIsLowX ? pointingCluster.GetInnerVertex().GetPosition() : pointingCluster.GetOuterVertex().GetPosition());
        highXEnd = (innerIsLowX ? pointingCluster.GetOuterVertex().GetPosition() : pointingCluster.GetInnerVertex().GetPosition());
    }
    catch (const StatusCodeException &)
    {
        return false;
    }

    const CartesianVector lowXUnitVector((lowXEnd - splitPosition).GetUnitVector());
    const CartesianVector highXUnitVector((highXEnd - splitPosition).GetUnitVector());

    CaloHitList caloHitList;
    pCurrentCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    std::string originalListName, fragmentListName;
    const ClusterList clusterList(1, pCurrentCluster);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, originalListName, fragmentListName));

    pLowXCluster = nullptr;
    pHighXCluster = nullptr;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector unitVector((pCaloHit->GetPositionVector() - splitPosition).GetUnitVector());
        const float dotProductLowX(unitVector.GetDotProduct(lowXUnitVector));
        const float dotProductHighX(unitVector.GetDotProduct(highXUnitVector));

        const Cluster *&pClusterToModify((dotProductLowX > dotProductHighX) ? pLowXCluster : pHighXCluster);

        if (!pClusterToModify)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
        }
    }

    if (!pLowXCluster || !pHighXCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        return false;
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool NViewTrackMatchingAlgorithm<T>::SortSplitPositions(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs)
{
    return (lhs.GetX() < rhs.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
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

    NViewMatchingAlgorithm<T>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    NViewMatchingAlgorithm<T>::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
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

template <typename T>
void NViewTrackMatchingAlgorithm<T>::PrepareInputClusters(ClusterList &preparedClusterList)
{
    for (ClusterList::iterator iter = preparedClusterList.begin(), iterEnd = preparedClusterList.end(); iter != iterEnd;)
    {
        try
        {
            this->AddToSlidingFitCache(*iter);
            ++iter;
        }
        catch (const StatusCodeException &statusCodeException)
        {
            preparedClusterList.erase(iter++);

            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::SetPfoParticleId(PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    pfoParameters.m_particleId = MU_MINUS; // Track
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::RemoveFromSlidingFitCache(const Cluster *const pCluster)
{
    TwoDSlidingFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NViewTrackMatchingAlgorithm<T>::TidyUp()
{
    m_slidingFitResultMap.clear();
    return NViewMatchingAlgorithm<T>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode NViewTrackMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    return NViewMatchingAlgorithm<T>::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class NViewTrackMatchingAlgorithm<TwoViewMatchingControl<float>>;
template class NViewTrackMatchingAlgorithm<TwoViewMatchingControl<TwoViewTransverseOverlapResult>>;
template class NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<float>>;
template class NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<TransverseOverlapResult>>;
template class NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<LongitudinalOverlapResult>>;
template class NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<FragmentOverlapResult>>;

} // namespace lar_content
