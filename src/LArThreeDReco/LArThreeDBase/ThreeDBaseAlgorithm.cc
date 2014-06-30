/**
 *  @file   LArContent/src/LArThreeDReco/LArThreeDBase/ThreeDBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimension algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"
#include "LArObjects/LArPointingCluster.h"
#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

template <typename T>
ThreeDBaseAlgorithm<T>::ThreeDBaseAlgorithm() :
    m_pInputClusterListU(NULL),
    m_pInputClusterListV(NULL),
    m_pInputClusterListW(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
ThreeDBaseAlgorithm<T>::~ThreeDBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool ThreeDBaseAlgorithm<T>::CreateThreeDParticles(const ProtoParticleVector &protoParticleVector)
{
    bool particlesMade(false);
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (typename ProtoParticleVector::const_iterator iter = protoParticleVector.begin(), iterEnd = protoParticleVector.end(); iter != iterEnd; ++iter)
    {
        // TODO - correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
        pfoParameters.m_clusterList.insert(iter->m_clusterListU.begin(), iter->m_clusterListU.end());
        pfoParameters.m_clusterList.insert(iter->m_clusterListV.begin(), iter->m_clusterListV.end());
        pfoParameters.m_clusterList.insert(iter->m_clusterListW.begin(), iter->m_clusterListW.end());

        ParticleFlowObject *pPfo(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
        particlesMade = true;
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool ThreeDBaseAlgorithm<T>::MakeClusterMerges(const ClusterMergeMap &clusterMergeMap)
{
    ClusterList deletedClusters;

    for (ClusterMergeMap::const_iterator iter = clusterMergeMap.begin(), iterEnd = clusterMergeMap.end(); iter != iterEnd; ++iter)
    {
        Cluster *pParentCluster = iter->first;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pParentCluster));
        const std::string clusterListName((TPC_VIEW_U == hitType) ? this->GetClusterListNameU() : (TPC_VIEW_V == hitType) ? this->GetClusterListNameV() : this->GetClusterListNameW());

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        for (ClusterList::const_iterator dIter = iter->second.begin(), dIterEnd = iter->second.end(); dIter != dIterEnd; ++dIter)
        {
            Cluster *pDaughterCluster = *dIter;

            if (deletedClusters.count(pParentCluster) || deletedClusters.count(pDaughterCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster, clusterListName, clusterListName));
            this->UpdateUponMerge(pParentCluster, pDaughterCluster);
            deletedClusters.insert(pDaughterCluster);
        }
    }

    return !(deletedClusters.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool ThreeDBaseAlgorithm<T>::MakeClusterSplits(const SplitPositionMap &splitPositionMap)
{
    bool changesMade(false);

    for (SplitPositionMap::const_iterator iter = splitPositionMap.begin(), iterEnd = splitPositionMap.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCurrentCluster = iter->first;
        CartesianPointList splitPositions(iter->second);
        std::sort(splitPositions.begin(), splitPositions.end(), ThreeDBaseAlgorithm::SortSplitPositions);

        const HitType hitType(LArClusterHelper::GetClusterHitType(pCurrentCluster));
        const std::string clusterListName((TPC_VIEW_U == hitType) ? this->GetClusterListNameU() : (TPC_VIEW_V == hitType) ? this->GetClusterListNameV() : this->GetClusterListNameW());

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        for (CartesianPointList::const_iterator sIter = splitPositions.begin(), sIterEnd = splitPositions.end(); sIter != sIterEnd; ++sIter)
        {
            Cluster *pLowXCluster(NULL), *pHighXCluster(NULL);
            this->MakeClusterSplit(*sIter, pCurrentCluster, pLowXCluster, pHighXCluster);

            this->UpdateUponSplit(pLowXCluster, pHighXCluster, pCurrentCluster);
            changesMade = true;
            pCurrentCluster = pHighXCluster;
        }
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::MakeClusterSplit(const CartesianVector &splitPosition, Cluster *&pCurrentCluster, Cluster *&pLowXCluster, Cluster *&pHighXCluster) const
{
    pLowXCluster = NULL;
    pHighXCluster = NULL;

    std::string originalListName, fragmentListName;
    ClusterList clusterList; clusterList.insert(pCurrentCluster);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, originalListName, fragmentListName));

    CaloHitList caloHitList;
    pCurrentCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);

    LArPointingCluster pointingCluster(pCurrentCluster);
    const bool innerIsLowX(pointingCluster.GetInnerVertex().GetPosition().GetX() < pointingCluster.GetOuterVertex().GetPosition().GetX());
    const CartesianVector &lowXEnd(innerIsLowX ? pointingCluster.GetInnerVertex().GetPosition() : pointingCluster.GetOuterVertex().GetPosition());
    const CartesianVector &highXEnd(innerIsLowX ? pointingCluster.GetOuterVertex().GetPosition() : pointingCluster.GetInnerVertex().GetPosition());

    const CartesianVector lowXUnitVector((lowXEnd -splitPosition).GetUnitVector());
    const CartesianVector highXUnitVector((highXEnd -splitPosition).GetUnitVector());

    for (CaloHitList::const_iterator hIter = caloHitList.begin(), hIterEnd = caloHitList.end(); hIter != hIterEnd; ++hIter)
    {
        CaloHit *pCaloHit = *hIter;
        const CartesianVector unitVector((pCaloHit->GetPositionVector() - splitPosition).GetUnitVector());

        const float dotProductLowX(unitVector.GetDotProduct(lowXUnitVector));
        const float dotProductHighX(unitVector.GetDotProduct(highXUnitVector));
        Cluster *&pClusterToModify((dotProductLowX > dotProductHighX) ? pLowXCluster : pHighXCluster);

        if (NULL == pClusterToModify)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.insert(pCaloHit);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
        }
    }

    if ((NULL == pLowXCluster) || (NULL == pHighXCluster))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool ThreeDBaseAlgorithm<T>::SortSplitPositions(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs)
{
    return (lhs.GetX() < rhs.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::UpdateUponMerge(Cluster *const pEnlargedCluster, Cluster *const pDeletedCluster)
{
    this->UpdateUponDeletion(pDeletedCluster);
    this->UpdateUponDeletion(pEnlargedCluster);
    this->UpdateForNewCluster(pEnlargedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::UpdateUponSplit(Cluster *const pSplitCluster1, Cluster *const pSplitCluster2, Cluster *const pDeletedCluster)
{
    this->UpdateUponDeletion(pDeletedCluster);
    this->UpdateForNewCluster(pSplitCluster1);
    this->UpdateForNewCluster(pSplitCluster2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::UpdateForNewCluster(Cluster *const pNewCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pNewCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((TPC_VIEW_U == hitType) ? m_clusterListU : (TPC_VIEW_V == hitType) ? m_clusterListV : m_clusterListW);

    if (!clusterList.insert(pNewCluster).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    const ClusterList &clusterList1((TPC_VIEW_U == hitType) ? m_clusterListV : m_clusterListU);
    const ClusterList &clusterList2((TPC_VIEW_W == hitType) ? m_clusterListV : m_clusterListW);

    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iter1End = clusterList1.end(); iter1 != iter1End; ++iter1)
    {
        for (ClusterList::const_iterator iter2 = clusterList2.begin(), iter2End = clusterList2.end(); iter2 != iter2End; ++iter2)
        {
            if (TPC_VIEW_U == hitType)
            {
                this->CalculateOverlapResult(pNewCluster, *iter1, *iter2);
            }
            else if (TPC_VIEW_V == hitType)
            {
                this->CalculateOverlapResult(*iter1, pNewCluster, *iter2);
            }
            else
            {
                this->CalculateOverlapResult(*iter1, *iter2, pNewCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::UpdateUponDeletion(Cluster *const pDeletedCluster)
{
    ClusterList::iterator iterU = m_clusterListU.find(pDeletedCluster);
    ClusterList::iterator iterV = m_clusterListV.find(pDeletedCluster);
    ClusterList::iterator iterW = m_clusterListW.find(pDeletedCluster);

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
void ThreeDBaseAlgorithm<T>::RemoveUnavailableTensorElements()
{
    const typename TensorType::ClusterNavigationMap &navigationMapUV(m_overlapTensor.GetClusterNavigationMapUV());
    const typename TensorType::ClusterNavigationMap &navigationMapVW(m_overlapTensor.GetClusterNavigationMapVW());
    const typename TensorType::ClusterNavigationMap &navigationMapWU(m_overlapTensor.GetClusterNavigationMapWU());
    ClusterList usedClusters;

    for (typename TensorType::ClusterNavigationMap::const_iterator iter = navigationMapUV.begin(), iterEnd = navigationMapUV.end(); iter != iterEnd; ++iter)
    {
        if (!(iter->first->IsAvailable()))
            usedClusters.insert(iter->first);
    }

    for (typename TensorType::ClusterNavigationMap::const_iterator iter = navigationMapVW.begin(), iterEnd = navigationMapVW.end(); iter != iterEnd; ++iter)
    {
        if (!(iter->first->IsAvailable()))
            usedClusters.insert(iter->first);
    }

    for (typename TensorType::ClusterNavigationMap::const_iterator iter = navigationMapWU.begin(), iterEnd = navigationMapWU.end(); iter != iterEnd; ++iter)
    {
        if (!(iter->first->IsAvailable()))
            usedClusters.insert(iter->first);
    }

    for (ClusterList::const_iterator iter = usedClusters.begin(), iterEnd = usedClusters.end(); iter != iterEnd; ++iter)
    {
        this->UpdateUponDeletion(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::SelectAllInputClusters()
{
    this->SelectInputClusters(m_pInputClusterListU, m_clusterListU);
    this->SelectInputClusters(m_pInputClusterListV, m_clusterListV);
    this->SelectInputClusters(m_pInputClusterListW, m_clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::PreparationStep()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::TidyUp()
{
    m_overlapTensor.Clear();

    m_pInputClusterListU = NULL;
    m_pInputClusterListV = NULL;
    m_pInputClusterListW = NULL;

    m_clusterListU.clear();
    m_clusterListV.clear();
    m_clusterListW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode ThreeDBaseAlgorithm<T>::Run()
{
    try
    {
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
            m_inputClusterListNameU, m_pInputClusterListU));
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
            m_inputClusterListNameV, m_pInputClusterListV));
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
            m_inputClusterListNameW, m_pInputClusterListW));

        if ((NULL == m_pInputClusterListU) || (NULL == m_pInputClusterListV) || (NULL == m_pInputClusterListW))
        {
            std::cout << "ThreeDBaseAlgorithm: one or more input cluster lists unavailable." << std::endl;
            throw StatusCodeException(STATUS_CODE_SUCCESS);
        }

        this->SelectAllInputClusters();
        this->PreparationStep();
        this->PerformMainLoop();
        this->ExamineTensor();
        this->TidyUp();
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();

        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::PerformMainLoop()
{
    for (ClusterList::const_iterator iterU = m_clusterListU.begin(), iterUEnd = m_clusterListU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterList::const_iterator iterV = m_clusterListV.begin(), iterVEnd = m_clusterListV.end(); iterV != iterVEnd; ++iterV)
        {
            for (ClusterList::const_iterator iterW = m_clusterListW.begin(), iterWEnd = m_clusterListW.end(); iterW != iterWEnd; ++iterW)
                this->CalculateOverlapResult(*iterU, *iterV, *iterW);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode ThreeDBaseAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

template class ThreeDBaseAlgorithm<float>;
template class ThreeDBaseAlgorithm<TransverseOverlapResult>;
template class ThreeDBaseAlgorithm<LongitudinalOverlapResult>;
template class ThreeDBaseAlgorithm<FragmentOverlapResult>;

} // namespace lar
