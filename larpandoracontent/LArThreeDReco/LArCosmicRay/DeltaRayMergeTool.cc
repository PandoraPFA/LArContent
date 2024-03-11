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
    m_maxDRSeparationFromTrack(1.5f),
    m_maxClusterSeparation(3.f),
    m_maxVertexSeparation(10.f),
    m_maxGoodMatchReducedChiSquared(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    m_pParentAlgorithm = pAlgorithm;

    if (PandoraContentApi::GetSettings(*m_pParentAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    return this->ExamineConnectedElements(overlapTensor);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::ExamineConnectedElements(TensorType &overlapTensor) const
{
    bool mergeMade(false), mergesMade(false), finishedTwoViewMerges(false);

    do
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
                usedKeyClusters.insert(element.GetClusterU());

            if (elementList.size() < 2)
                continue;

            if (!finishedTwoViewMerges && this->MakeTwoCommonViewMerges(elementList))
            {
                mergeMade = true;
                mergesMade = true;
                break;
            }

            finishedTwoViewMerges = true;

            if (this->MakeOneCommonViewMerges(elementList))
            {
                mergeMade = true;
                mergesMade = true;
                break;
            }
        }
    } while (mergeMade);

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::MakeTwoCommonViewMerges(const TensorType::ElementList &elementList) const
{
    for (auto iter1 = elementList.begin(); iter1 != elementList.end(); ++iter1)
    {
        const TensorType::Element &element1(*iter1);

        for (auto iter2 = std::next(iter1); iter2 != elementList.end(); ++iter2)
        {
            const TensorType::Element &element2(*iter2);

            for (const HitType hitType1 : {TPC_VIEW_U, TPC_VIEW_V})
            {
                if ((element1.GetCluster(hitType1) == element2.GetCluster(hitType1)))
                {
                    for (const HitType hitType2 : {TPC_VIEW_V, TPC_VIEW_W})
                    {
                        if (hitType1 == hitType2)
                            continue;

                        if ((element1.GetCluster(hitType2) == element2.GetCluster(hitType2)))
                        {
                            const HitType mergeHitType(hitType1 == TPC_VIEW_U ? (hitType2 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V) : TPC_VIEW_U);
                            const Cluster *pClusterToEnlarge(element1.GetCluster(mergeHitType)), *pClusterToDelete(element2.GetCluster(mergeHitType));

                            if (this->AreAssociated(element1, element2, mergeHitType))
                            {
                                m_pParentAlgorithm->UpdateUponDeletion(pClusterToEnlarge);
                                m_pParentAlgorithm->UpdateUponDeletion(pClusterToDelete);

                                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                    PandoraContentApi::ReplaceCurrentList<Cluster>(
                                        *m_pParentAlgorithm, m_pParentAlgorithm->GetClusterListName(mergeHitType)));

                                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                    PandoraContentApi::MergeAndDeleteClusters(*m_pParentAlgorithm, pClusterToEnlarge, pClusterToDelete));

                                m_pParentAlgorithm->UpdateForNewClusters({pClusterToEnlarge}, {nullptr});

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
    const PfoList &commonMuonPfoList1(element1.GetOverlapResult().GetCommonMuonPfoList());
    const PfoList &commonMuonPfoList2(element2.GetOverlapResult().GetCommonMuonPfoList());

    // Demand the elements to have a shared common muon
    PfoList commonMuonPfoList;
    this->CombineCommonMuonPfoLists(commonMuonPfoList1, commonMuonPfoList2, commonMuonPfoList);

    if (commonMuonPfoList.empty())
        return false;

    const Cluster *const pCluster1(element1.GetCluster(mergeHitType)), *const pCluster2(element2.GetCluster(mergeHitType));

    PfoList connectedMuonPfoList1, connectedMuonPfoList2;
    this->GetConnectedMuons(pCluster1, commonMuonPfoList1, connectedMuonPfoList1);
    this->GetConnectedMuons(pCluster2, commonMuonPfoList2, connectedMuonPfoList2);

    if (connectedMuonPfoList1.empty() || connectedMuonPfoList2.empty())
        return (this->IsBrokenCluster(pCluster1, pCluster2));

    for (const ParticleFlowObject *const pConnectedMuon1 : connectedMuonPfoList1)
    {
        for (const ParticleFlowObject *const pConnectedMuon2 : connectedMuonPfoList2)
        {
            if ((pConnectedMuon1 == pConnectedMuon2) && (this->IsHiddenByTrack(pConnectedMuon1, pCluster1, pCluster2)))
                return true;
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

void DeltaRayMergeTool::GetConnectedMuons(const Cluster *const pDeltaRayCluster, const PfoList &commonMuonPfoList, PfoList &connectedMuonPfoList) const
{
    for (const ParticleFlowObject *const pCommonMuonPfo : commonMuonPfoList)
    {
        if (this->IsConnected(pDeltaRayCluster, pCommonMuonPfo))
            connectedMuonPfoList.push_back(pCommonMuonPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::IsConnected(const Cluster *const pCluster, const Pfo *const pCommonMuonPfo) const
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

bool DeltaRayMergeTool::IsHiddenByTrack(const ParticleFlowObject *const pMuonPfo, const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    CaloHitList vertices1, vertices2;
    this->FindVertices(pMuonPfo, pCluster1, vertices1);
    this->FindVertices(pMuonPfo, pCluster2, vertices2);

    if (vertices1.empty() || vertices2.empty())
        return false;

    for (const CaloHit *const pCaloHit : vertices1)
    {
        const float separation(LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), vertices2));

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

bool DeltaRayMergeTool::MakeOneCommonViewMerges(const TensorType::ElementList &elementList) const
{
    for (auto iter1 = elementList.begin(); iter1 != elementList.end(); ++iter1)
    {
        const TensorType::Element &element1(*iter1);

        for (auto iter2 = std::next(iter1); iter2 != elementList.end(); ++iter2)
        {
            const TensorType::Element &element2(*iter2);

            for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                if (element1.GetCluster(hitType) == element2.GetCluster(hitType))
                {
                    const HitType mergeHitType1(hitType == TPC_VIEW_U ? TPC_VIEW_V : hitType == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                    const HitType mergeHitType2(mergeHitType1 == TPC_VIEW_U ? TPC_VIEW_V
                            : mergeHitType1 == TPC_VIEW_V                   ? TPC_VIEW_W
                                                                            : TPC_VIEW_U);

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

                    float reducedChiSquared(std::numeric_limits<float>::max());
                    StatusCode status(m_pParentAlgorithm->PerformThreeViewMatching(caloHitList1, caloHitList2, caloHitList3, reducedChiSquared));

                    if (status == STATUS_CODE_NOT_FOUND)
                        continue;

                    if (reducedChiSquared < m_maxGoodMatchReducedChiSquared)
                    {
                        m_pParentAlgorithm->UpdateUponDeletion(pClusterToEnlarge1);
                        m_pParentAlgorithm->UpdateUponDeletion(pClusterToDelete1);

                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraContentApi::ReplaceCurrentList<Cluster>(*m_pParentAlgorithm, m_pParentAlgorithm->GetClusterListName(mergeHitType1)));

                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraContentApi::MergeAndDeleteClusters(*m_pParentAlgorithm, pClusterToEnlarge1, pClusterToDelete1));

                        m_pParentAlgorithm->UpdateForNewClusters({pClusterToEnlarge1}, {nullptr});

                        m_pParentAlgorithm->UpdateUponDeletion(pClusterToEnlarge2);
                        m_pParentAlgorithm->UpdateUponDeletion(pClusterToDelete2);

                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraContentApi::ReplaceCurrentList<Cluster>(*m_pParentAlgorithm, m_pParentAlgorithm->GetClusterListName(mergeHitType2)));

                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraContentApi::MergeAndDeleteClusters(*m_pParentAlgorithm, pClusterToEnlarge2, pClusterToDelete2));

                        m_pParentAlgorithm->UpdateForNewClusters({pClusterToEnlarge2}, {nullptr});

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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDRSeparationFromTrack", m_maxDRSeparationFromTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxVertexSeparation", m_maxVertexSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxGoodMatchReducedChiSquared", m_maxGoodMatchReducedChiSquared));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
