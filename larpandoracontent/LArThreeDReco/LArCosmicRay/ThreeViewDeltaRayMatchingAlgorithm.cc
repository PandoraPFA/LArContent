/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

ThreeViewDeltaRayMatchingAlgorithm::ThreeViewDeltaRayMatchingAlgorithm() :
    m_nMaxTensorToolRepeats(10),
    m_minClusterCaloHits(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::GetConnectedElements(const Cluster *const pCluster1, const bool hasAssociatedMuon, TensorType::ElementList &elementList, ClusterSet &checkedClusters)
{
    if (checkedClusters.count(pCluster1))
        return;

    if (!pCluster1->IsAvailable())
        return;

    auto &theTensor(this->GetMatchingControl().GetOverlapTensor());
    
    const HitType &hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    
    if (hitType1 == TPC_VIEW_U)
        checkedClusters.insert(pCluster1);
    
    auto &navigationMap1((hitType1 == TPC_VIEW_U) ? theTensor.GetClusterNavigationMapUV() : (hitType1 == TPC_VIEW_V) ? theTensor.GetClusterNavigationMapVW() :
        theTensor.GetClusterNavigationMapWU());

    auto iter1 = navigationMap1.find(pCluster1);

    if (iter1 == navigationMap1.end())
        throw StatusCodeException(STATUS_CODE_FAILURE);
    
    for (const Cluster *const pCluster2 : iter1->second)
    {
        if (checkedClusters.count(pCluster2))
            continue;

        if (!pCluster2->IsAvailable())
            continue;

        const HitType &hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
        auto &navigationMap2((hitType2 == TPC_VIEW_U) ? theTensor.GetClusterNavigationMapUV() : (hitType2 == TPC_VIEW_V) ? theTensor.GetClusterNavigationMapVW() :
            theTensor.GetClusterNavigationMapWU());

        auto iter2 = navigationMap2.find(pCluster2);

        if (iter2 == navigationMap2.end())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        for (const Cluster *const pCluster3 : iter2->second)
        {
            if (checkedClusters.count(pCluster3))
                continue;

            if (!pCluster3->IsAvailable())
                continue;

            const Cluster *const pClusterU((hitType1 == TPC_VIEW_U) ? pCluster1 : (hitType2 == TPC_VIEW_U) ? pCluster2 : pCluster3);
            const Cluster *const pClusterV((hitType1 == TPC_VIEW_V) ? pCluster1 : (hitType2 == TPC_VIEW_V) ? pCluster2 : pCluster3);
            const Cluster *const pClusterW((hitType1 == TPC_VIEW_W) ? pCluster1 : (hitType2 == TPC_VIEW_W) ? pCluster2 : pCluster3);

            if ((pClusterU == pClusterV) || (pClusterV == pClusterW) || (pClusterU == pClusterW))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // now find in the tensor
            try
            {
                auto &overlapResult(theTensor.GetOverlapResult(pClusterU, pClusterV, pClusterW));
            
                const PfoList &commonMuonPfoList(overlapResult.GetCommonMuonPfoList());

                if (!hasAssociatedMuon && commonMuonPfoList.size())
                    continue;

                if (hasAssociatedMuon && commonMuonPfoList.empty())
                    continue;
                
                bool found = false;
                for (const TensorType::Element &t : elementList)
                {
                    if ((t.GetClusterU() == pClusterU) && (t.GetClusterV() == pClusterV) && (t.GetClusterW() == pClusterW))
                        found = true;
                }

                if (!found)
                {
                    TensorType::Element element(pClusterU, pClusterV, pClusterW, overlapResult);
                    elementList.push_back(element);
                }
                

                //TensorType::Element element(pClusterU, pClusterV, pClusterW, overlapResult);
                //elementList.push_back(element);

                // check if we really need both
                this->GetConnectedElements(pCluster2, hasAssociatedMuon, elementList, checkedClusters);
                this->GetConnectedElements(pCluster3, hasAssociatedMuon, elementList, checkedClusters);
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }
            
    }

    std::sort(elementList.begin(), elementList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void ThreeViewDeltaRayMatchingAlgorithm::GetUnambiguousElements(const bool hasAssociatedMuon, TensorType::ElementList &elementList)
{
    auto &theMatrix(this->GetMatchingControl().GetOverlapTensor());

    for (auto iterU = theMatrix.begin(); iterU != theMatrix.end(); ++iterU)
    {
        ClusterSet checkedClusters;
        TensorType::ElementList tempElementList;
        this->GetConnectedElements(iterU->first, hasAssociatedMuon, tempElementList, checkedClusters);

        if (tempElementList.size() != 1)
            continue;

        TensorType::Element &unambiguouslement(tempElementList.front());
        const Cluster *const pClusterU(unambiguouslement.GetClusterU()), *const pClusterV(unambiguouslement.GetClusterV()), *const pClusterW(unambiguouslement.GetClusterW());
        
        // ATTN With HIT_CUSTOM definitions, it is possible to navigate from different U clusters to same combination
        if (iterU->first != pClusterU)
            continue;

        if (!pClusterU || !pClusterV || !pClusterW)
            continue;

        auto iterV = iterU->second.find(pClusterV);
        if (iterU->second.end() == iterV)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        auto iterW = iterV->second.find(pClusterW);
        if (iterV->second.end() == iterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        TensorType::Element element(pClusterU, pClusterV, pClusterW, iterW->second);
        elementList.push_back(element);
    }

    std::sort(elementList.begin(), elementList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewDeltaRayMatchingAlgorithm::DoesClusterPassTesorThreshold(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
        return false;

    return true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    DeltaRayOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult));

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    DeltaRayOverlapResult &overlapResult) const
{
    PfoList commonMuonPfoList;
    this->FindCommonMuonParents(pClusterU, pClusterV, pClusterW, commonMuonPfoList);

    if (commonMuonPfoList.empty())
        return STATUS_CODE_NOT_FOUND;
    
    float chiSquaredSum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    XOverlap xOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    
    StatusCode statusCode(this->PerformThreeViewMatching(pClusterU, pClusterV, pClusterW, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xOverlapObject));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

     overlapResult = DeltaRayOverlapResult(nMatchedSamplingPoints, nSamplingPoints, chiSquaredSum, xOverlapObject, commonMuonPfoList);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::FindCommonMuonParents(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    PfoList &commonMuonPfoList) const
{
    ClusterList consideredClustersU, consideredClustersV, consideredClustersW;
    PfoList nearbyMuonPfosU, nearbyMuonPfosV, nearbyMuonPfosW;
    this->GetNearbyMuonPfos(pClusterU, consideredClustersU, nearbyMuonPfosU);
    this->GetNearbyMuonPfos(pClusterV, consideredClustersV, nearbyMuonPfosV);
    this->GetNearbyMuonPfos(pClusterW, consideredClustersW, nearbyMuonPfosW);

    for (const ParticleFlowObject *const pNearbyMuonU : nearbyMuonPfosU)
    {
        for (const ParticleFlowObject *const pNearbyMuonV : nearbyMuonPfosV)
        {
            for (const ParticleFlowObject *const pNearbyMuonW : nearbyMuonPfosW)
            {
                if ((pNearbyMuonU == pNearbyMuonV) && (pNearbyMuonV == pNearbyMuonW))
                {
                    commonMuonPfoList.push_back(pNearbyMuonU);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (TensorToolVector::const_iterator toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end(); )
    {
        if ((*toolIter)->Run(this, this->GetMatchingControl().GetOverlapTensor()))
        {
            toolIter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++toolIter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ClusterRebuilding", m_reclusteringAlgorithmName));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName));
    
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "DeltaRayTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DeltaRayTensorTool *const pDeltaRayTensorTool(dynamic_cast<DeltaRayTensorTool*>(*iter));

        if (!pDeltaRayTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayTensorTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));
    
    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content

