/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMergeTool.cc
 *
 *  @brief  Implementation of the delta ray merge tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMergeTool.h"

using namespace pandora;

namespace lar_content
{

DeltaRayMergeTool::DeltaRayMergeTool() :
    m_maxUnambiguousClusterSeparation(1.f),
    m_maxDRSeparationFromTrack(1.5f),
    m_maxVertexSeparation(10.f),
    m_maxClusterSeparation(3.f),
    m_maxGoodMatchReducedChiSquared(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool mergesMade(false);

    this->ExamineConnectedElements(pAlgorithm, overlapTensor, mergesMade);

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor, bool &mergesMade) const
{
    bool mergeMade(true);

    while (mergeMade)
    {
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            TensorType::ElementList elementList;
            overlapTensor.GetConnectedElements(pKeyCluster, true, elementList);
            
            for (const TensorType::Element &element : elementList)
            {
                if (usedKeyClusters.count(element.GetClusterU()))
                    continue;

                usedKeyClusters.insert(element.GetClusterU());
            }

            if (elementList.size() < 2)
                continue;

            if (this->MakeTwoCommonViewMerges(pAlgorithm, elementList))
            {
                mergeMade = true; mergesMade = true;
                break;
            }
        }
    }

    mergeMade = true;

    while (mergeMade)
    {
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            TensorType::ElementList elementList;
            overlapTensor.GetConnectedElements(pKeyCluster, true, elementList);

            for (const TensorType::Element &element : elementList)
            {
                if (usedKeyClusters.count(element.GetClusterU()))
                    continue;

                usedKeyClusters.insert(element.GetClusterU());
            }

            if (elementList.size() < 2)
                continue;

            if (this->MakeOneCommonViewMerges(pAlgorithm, elementList))
            {
                mergeMade = true; mergesMade = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::MakeTwoCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const
{
    const HitTypeVector hitTypeVector1({TPC_VIEW_U, TPC_VIEW_V});
    const HitTypeVector hitTypeVector2({TPC_VIEW_V, TPC_VIEW_W});

    for (const TensorType::Element &element1 : elementList)
    {
        for (const TensorType::Element &element2 : elementList)
        {
            if ((element1.GetCluster(TPC_VIEW_U) == element2.GetCluster(TPC_VIEW_U)) && (element1.GetCluster(TPC_VIEW_V) == element2.GetCluster(TPC_VIEW_V)) &&
                (element1.GetCluster(TPC_VIEW_W) == element2.GetCluster(TPC_VIEW_W)))
            {
                continue;
            }

            for (const HitType &hitType1 : hitTypeVector1)
            {
                if ((element1.GetCluster(hitType1) == element2.GetCluster(hitType1)))
                {
                    for(const HitType &hitType2 : hitTypeVector2)
                    {
                        if (hitType1 == hitType2)
                            continue;
                        
                        if ((element1.GetCluster(hitType2) == element2.GetCluster(hitType2)))
                        {
                            const HitType mergeHitType(hitType1 == TPC_VIEW_U ? (hitType2 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V) : TPC_VIEW_U);
                            const Cluster *pClusterToEnlarge(element1.GetCluster(mergeHitType)), *pClusterToDelete(element2.GetCluster(mergeHitType));

                            if (this->AreAssociated(element1, element2, mergeHitType))
                            {
                                pAlgorithm->UpdateUponDeletion(pClusterToEnlarge); pAlgorithm->UpdateUponDeletion(pClusterToDelete);
                                    
                                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                                   pAlgorithm->GetClusterListName(LArClusterHelper::GetClusterHitType(pClusterToEnlarge))));
                                
                                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge, pClusterToDelete));
                                    
                                pAlgorithm->UpdateForNewClusters({pClusterToEnlarge}, {nullptr});

                                return true; 
                            }
                        }
                    }
                }
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::AreAssociated(const TensorType::Element &element1, const TensorType::Element &element2, const HitType &mergeHitType) const
{
    // Demand the elements to have a shared common muon
    PfoList commonMuonPfoList;
    this->CombineCommonMuonPfoLists(element1.GetOverlapResult().GetCommonMuonPfoList(), element2.GetOverlapResult().GetCommonMuonPfoList(), commonMuonPfoList);

    if (commonMuonPfoList.empty())
        return false;
    
    const Cluster *const pCluster1(element1.GetCluster(mergeHitType)), *const pCluster2(element2.GetCluster(mergeHitType));

    PfoList connectedMuonPfoList1, connectedMuonPfoList2;
    this->GetConnectedMuons(element1.GetOverlapResult().GetCommonMuonPfoList(), pCluster1, connectedMuonPfoList1);
    this->GetConnectedMuons(element2.GetOverlapResult().GetCommonMuonPfoList(), pCluster2, connectedMuonPfoList2);

    if (connectedMuonPfoList1.empty() || connectedMuonPfoList2.empty())
    {
      if (this->IsBrokenCluster(pCluster1, pCluster2))                                                                                                                                   
          return true;
    }

    for (const ParticleFlowObject *const pConnectedMuon1 : connectedMuonPfoList1)
    {
        for (const ParticleFlowObject *const pConnectedMuon2 : connectedMuonPfoList2)
        {
            if (pConnectedMuon1 == pConnectedMuon2)
            {
                if (this->IsHiddenTrack(pConnectedMuon1, pCluster1, pCluster2))
                    return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::CombineCommonMuonPfoLists(const PfoList &commonMuonPfoList1, const PfoList &commonMuonPfoList2, PfoList &commonMuonPfoList) const
{
    for (const ParticleFlowObject *const pCommonMuonPfo1 : commonMuonPfoList1)
    {
        for (const ParticleFlowObject *const pCommonMuonPfo2 : commonMuonPfoList2)
        {
            if (pCommonMuonPfo1 == pCommonMuonPfo2)
                commonMuonPfoList.push_back(pCommonMuonPfo1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::GetConnectedMuons(const PfoList &commonMuonPfoList, const Cluster *const pClusterToEnlarge, PfoList &connectedMuonPfoList) const 
{
    for (const ParticleFlowObject *const pCommonMuonPfo : commonMuonPfoList)
    {
        if (this->IsConnected(pCommonMuonPfo, pClusterToEnlarge))
            connectedMuonPfoList.push_back(pCommonMuonPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool DeltaRayMergeTool::IsConnected(const Pfo *const pCommonMuonPfo, const Cluster *const pCluster) const
{
    HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pCommonMuonPfo, hitType, muonClusterList);    

    if (muonClusterList.size() != 1)
        return false;
    
    const float separation(LArClusterHelper::GetClosestDistance(pCluster, muonClusterList));

    return separation < m_maxDRSeparationFromTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::IsBrokenCluster(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete) const
{
    const float clusterSeparation(LArClusterHelper::GetClosestDistance(pClusterToEnlarge, pClusterToDelete));

    return clusterSeparation < m_maxClusterSeparation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::IsHiddenTrack(const ParticleFlowObject *const pMuonPfo, const Cluster *const pCluster1, const Cluster *const pCluster2) const 
{
    CaloHitList vertices1, vertices2;
    this->FindVertices(pMuonPfo, pCluster1, vertices1);
    this->FindVertices(pMuonPfo, pCluster2, vertices2); 

    for (const CaloHit *const pCaloHit : vertices1)
    {
        const float separation(LArMuonLeadingHelper::GetClosestDistance(pCaloHit, vertices2));
        
        if (separation < m_maxVertexSeparation)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::FindVertices(const Pfo *const pCommonMuonPfo, const Cluster *const pCluster, CaloHitList &vertexList) const
{
    HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pCommonMuonPfo, hitType, muonClusterList);    

    if (muonClusterList.size() != 1)
        return;
    
    CaloHitList caloHitList;
    muonClusterList.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pCluster) < m_maxDRSeparationFromTrack)
            vertexList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::MakeOneCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    for (const TensorType::Element &element1 : elementList)
    {
        for (const TensorType::Element &element2 : elementList)
        {
            if ((element1.GetCluster(TPC_VIEW_U) == element2.GetCluster(TPC_VIEW_U)) && (element1.GetCluster(TPC_VIEW_V) == element2.GetCluster(TPC_VIEW_V)) &&
                (element1.GetCluster(TPC_VIEW_W) == element2.GetCluster(TPC_VIEW_W)))
            {
                continue;
            }

            for (const HitType &hitType : hitTypeVector)
            {
                if (element1.GetCluster(hitType) == element2.GetCluster(hitType))
                {
                    const HitType mergeHitType1(hitType == TPC_VIEW_U ? TPC_VIEW_V : hitType == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                    const HitType mergeHitType2(mergeHitType1 == TPC_VIEW_U ? TPC_VIEW_V : mergeHitType1 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                        
                    const Cluster *pClusterToEnlarge1 = element1.GetCluster(mergeHitType1), *pClusterToDelete1 = element2.GetCluster(mergeHitType1);
                    const Cluster *pClusterToEnlarge2 = element1.GetCluster(mergeHitType2), *pClusterToDelete2 = element2.GetCluster(mergeHitType2);

                    if ((pClusterToEnlarge1 == pClusterToDelete1) || (pClusterToEnlarge2 == pClusterToDelete2))
                        continue;
                        
                    if (!this->AreAssociated(element1, element2, mergeHitType1))
                        continue;

                    if (!this->AreAssociated(element1, element2, mergeHitType2))
                        continue;

                    CaloHitList caloHitList1, caloHitList2, caloHitList3;
                    pClusterToEnlarge1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
                    pClusterToDelete1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
                    pClusterToEnlarge2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
                    pClusterToDelete2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
                    element1.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);

                    XOverlap xOverlapObject(0.f,0.f,0.f,0.f,0.f, 0.f,0.f);
                    float chiSquaredSum(0.f);
                    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
                        
                    StatusCode status(pAlgorithm->PerformThreeViewMatching(caloHitList1, caloHitList2, caloHitList3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xOverlapObject));
                        
                    if (status == STATUS_CODE_NOT_FOUND)
                        continue;

                    const float reducedChiSquared(chiSquaredSum / nSamplingPoints);
                        
                    if (reducedChiSquared < m_maxGoodMatchReducedChiSquared)
                    {
                        pAlgorithm->UpdateUponDeletion(pClusterToEnlarge1); pAlgorithm->UpdateUponDeletion(pClusterToDelete1);

                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                            pAlgorithm->GetClusterListName(mergeHitType1)));
                        
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge1, pClusterToDelete1));

                        pAlgorithm->UpdateForNewClusters({pClusterToEnlarge1}, {nullptr});

                        pAlgorithm->UpdateUponDeletion(pClusterToEnlarge2); pAlgorithm->UpdateUponDeletion(pClusterToDelete2);

                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                            pAlgorithm->GetClusterListName(mergeHitType2)));
                        
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge2, pClusterToDelete2));
			   
                        pAlgorithm->UpdateForNewClusters({pClusterToEnlarge2}, {nullptr});

                        return true;
                    }
                }
            }
        }
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMergeTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDRSeparationFromTrack", m_maxDRSeparationFromTrack));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexSeparation", m_maxVertexSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterSeparation", m_maxClusterSeparation));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxUnambiguousClusterSeparation", m_maxUnambiguousClusterSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGoodMatchReducedChiSquared", m_maxGoodMatchReducedChiSquared));     
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
