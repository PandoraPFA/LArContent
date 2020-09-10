/**
 *  @file   larpandoracontent/LArVertex/DirectionClusterSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArDirection/DirectionClusterSplittingAlgorithm.h"

#include <utility>
#include <limits>

using namespace pandora;

namespace lar_content
{

DirectionClusterSplittingAlgorithm::DirectionClusterSplittingAlgorithm() :
    m_minClusterCaloHits(50),
    m_minClusterLengthSquared(30.f * 30.f),
    m_enableDirection(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionClusterSplittingAlgorithm::Run()
{
    if (!m_enableDirection)
        return STATUS_CODE_SUCCESS;

    //this->CountClusters();

    try
    {
        ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;
        this->SelectClusters(clusterVectorU, clusterVectorV, clusterVectorW);

        std::string originalListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalListName));

        for (auto pCluster : clusterVectorW)
        {
            try
            {
                TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetClusterDirection(pCluster);
                //fitResult.DrawFit();

                if (!fitResult.GetSplitObject().GetSplitApplied())
                    continue;

                this->SplitCluster(this->RetrieveSplitCaloHitClusterPair(pCluster, fitResult), TPC_VIEW_W); 
                this->SplitCluster(this->FindMatchingCaloHitClusterPair(fitResult, clusterVectorV), TPC_VIEW_V); 
                this->SplitCluster(this->FindMatchingCaloHitClusterPair(fitResult, clusterVectorU), TPC_VIEW_U); 
            }
            catch (...)
            {
                continue;
            }
            
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalListName));
    }
    catch (StatusCodeException &statusCodeException)
    {
        throw statusCodeException;
    }

    //this->CountClusters();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionClusterSplittingAlgorithm::CountClusters()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if (hitType == TPC_VIEW_U)
            std::cout << "Number U clusters: " << pClusterList->size() << std::endl;
        if (hitType == TPC_VIEW_V)
            std::cout << "Number V clusters: " << pClusterList->size() << std::endl;
        if (hitType == TPC_VIEW_W)
            std::cout << "Number W clusters: " << pClusterList->size() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionClusterSplittingAlgorithm::SelectClusters(ClusterVector &clusterVectorU, ClusterVector &clusterVectorV, ClusterVector &clusterVectorW)
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "DirectionClusterSplittingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterVector &selectedClusterVector((TPC_VIEW_U == hitType) ? clusterVectorU : (TPC_VIEW_V == hitType) ? clusterVectorV : clusterVectorW);

        if (!selectedClusterVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
        std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pCluster : sortedClusters)
        {   
            if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
                continue;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
                continue;

            selectedClusterVector.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Cluster* DirectionClusterSplittingAlgorithm::GetTargetCluster(ClusterVector &clusterVectorW)
{
    float currentLength(std::numeric_limits<float>::min());
    pandora::ClusterVector longestClusterList;

    for (const Cluster *const pCluster : clusterVectorW)
    {    
        std::cout << currentLength << std::endl;
        if (LArClusterHelper::GetLength(pCluster) > currentLength)
        {    
            currentLength = LArClusterHelper::GetLength(pCluster);
            longestClusterList.clear();
            longestClusterList.push_back(pCluster);
        }    
    }    
    
    if (longestClusterList.size() == 0)
        return NULL;
    else
        return (*(longestClusterList.begin()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<const pandora::CaloHit*, const pandora::Cluster*> DirectionClusterSplittingAlgorithm::RetrieveSplitCaloHitClusterPair(const pandora::Cluster* pCluster, TrackDirectionTool::DirectionFitObject &fitResult)
{
    float splitPosition(fitResult.GetSplitObject().GetSplitPosition());

    float closestDistance(std::numeric_limits<float>::max());
    pandora::CaloHitVector targetCaloHitVector;

    for (auto &hitCharge : fitResult.GetHitChargeVector())
    {
        if (std::abs(hitCharge.GetLongitudinalPosition() - splitPosition) < closestDistance)
        {
            closestDistance = std::abs(hitCharge.GetLongitudinalPosition() - splitPosition);
            targetCaloHitVector.clear();
            targetCaloHitVector.push_back(hitCharge.GetCaloHit());
        }
    }

    return std::make_pair(targetCaloHitVector.front(), pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<const pandora::CaloHit*, const pandora::Cluster*> DirectionClusterSplittingAlgorithm::FindMatchingCaloHitClusterPair(TrackDirectionTool::DirectionFitObject &fitResult, ClusterVector &clusterVector)
{
    const pandora::CaloHit* leftHit((fitResult.GetBackwardsFitCharges().back()).GetCaloHit());
    const pandora::CaloHit* rightHit((fitResult.GetForwardsFitCharges().front()).GetCaloHit());
    
    const float targetX((leftHit->GetPositionVector().GetX() + rightHit->GetPositionVector().GetX())/2);
    //const float lowerZ(leftHit->GetPositionVector().GetZ() < rightHit->GetPositionVector().GetZ() ? leftHit->GetPositionVector().GetZ() : rightHit->GetPositionVector().GetZ());
    //const float upperZ(leftHit->GetPositionVector().GetZ() > rightHit->GetPositionVector().GetZ() ? leftHit->GetPositionVector().GetZ() : rightHit->GetPositionVector().GetZ());

    float closestDistance(std::numeric_limits<float>::max());
    pandora::CaloHitVector targetCaloHitVector;
    pandora::ClusterVector targetClusterVector;

    for (auto pCluster : clusterVector)
    {
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.FillCaloHitList(caloHitList);    

        for (auto pCaloHit : caloHitList)
        {    
            if (std::abs(pCaloHit->GetPositionVector().GetX() - targetX) < closestDistance) 
            {    
                closestDistance = std::abs(pCaloHit->GetPositionVector().GetX() - targetX); 

                targetCaloHitVector.clear();
                targetCaloHitVector.push_back(pCaloHit);

                targetClusterVector.clear();
                targetClusterVector.push_back(pCluster);
            }    
        }  
    }
    
    if (targetCaloHitVector.size() == 0 || targetClusterVector.size() == 0)
        throw STATUS_CODE_NOT_FOUND;

    return std::make_pair(targetCaloHitVector.front(), targetClusterVector.front());
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DirectionClusterSplittingAlgorithm::SplitCluster(std::pair<const pandora::CaloHit*, const pandora::Cluster*> hitClusterPair, pandora::HitType listHitType)
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if (hitType == listHitType)
        {
            const StatusCode listChangeStatusCode(PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));
            if (listChangeStatusCode != STATUS_CODE_SUCCESS)
                throw listChangeStatusCode;
        }
    }

    const pandora::CaloHit* pCaloHit(hitClusterPair.first);
    const pandora::Cluster* pCluster(hitClusterPair.second);

    // Split cluster into two CaloHit lists
    PandoraContentApi::Cluster::Parameters firstParameters, secondParameters;

    this->DivideCaloHits(pCluster, pCaloHit, firstParameters.m_caloHitList, secondParameters.m_caloHitList);

    if (firstParameters.m_caloHitList.empty() || secondParameters.m_caloHitList.empty())
        return STATUS_CODE_NOT_ALLOWED;

    // Begin cluster fragmentation operations
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    const Cluster *pFirstCluster(NULL), *pSecondCluster(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, firstParameters, pFirstCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, secondParameters, pSecondCluster));

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionClusterSplittingAlgorithm::DivideCaloHits(const pandora::Cluster* pCluster, const pandora::CaloHit* pTargetCaloHit, CaloHitList &caloHitList1, CaloHitList &caloHitList2)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);    

    pandora::CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());

    //Sort by Z with lambda function
    std::sort(caloHitVector.begin(), caloHitVector.end(), [](auto pCaloHit1, auto pCaloHit2)  {return pCaloHit1->GetPositionVector().GetZ() < pCaloHit2->GetPositionVector().GetZ();});

    auto pLowZCaloHit(caloHitVector.front());

    for (auto pCaloHit : caloHitVector) 
    {
       if ((pLowZCaloHit->GetPositionVector() - pCaloHit->GetPositionVector()).GetMagnitude() < (pLowZCaloHit->GetPositionVector() - pTargetCaloHit->GetPositionVector()).GetMagnitude()) 
           caloHitList1.push_back(pCaloHit);
        else    
           caloHitList2.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableDirection", m_enableDirection));

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
