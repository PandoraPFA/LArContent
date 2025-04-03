/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the two view delta ray matching class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoViewDeltaRayMatchingAlgorithm::TwoViewDeltaRayMatchingAlgorithm() :
    m_nMaxMatrixToolRepeats(10),
    m_minClusterCaloHits(3),
    m_maxDistanceFromPrediction(2.f),
    m_maxGoodMatchReducedChiSquared(1.f),
    m_minDistanceFromMuon(1.f),
    m_maxDistanceToCollected(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewDeltaRayMatchingAlgorithm::HitTypeVector TwoViewDeltaRayMatchingAlgorithm::GetHitTypeVector()
{
    HitTypeVector hitTypeVector;

    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const unsigned int hitTypeIndex(this->GetMatchingControl().GetHitTypeIndex(hitType));

        if ((hitTypeIndex != 1) && (hitTypeIndex != 2))
            continue;

        hitTypeVector.push_back(hitType);
    }

    return hitTypeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TwoViewDeltaRayMatchingAlgorithm::GetCluster(const MatrixType::Element &element, const HitType hitType)
{
    const unsigned int hitTypeIndex(this->GetMatchingControl().GetHitTypeIndex(hitType));

    if ((hitTypeIndex != 1) && (hitTypeIndex != 2))
        return element.GetOverlapResult().GetBestMatchedAvailableCluster();

    return hitTypeIndex == 1 ? element.GetCluster1() : element.GetCluster2();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMatchingAlgorithm::DoesClusterPassTensorThreshold(const Cluster *const pCluster) const
{
    return (pCluster->GetNCaloHits() >= m_minClusterCaloHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const)
{
    TwoViewDeltaRayOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pCluster1, pCluster2, overlapResult));

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(
    const Cluster *const pCluster1, const Cluster *const pCluster2, TwoViewDeltaRayOverlapResult &overlapResult) const
{
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    PfoList commonMuonPfoList;
    this->FindCommonMuonParents(pCluster1, pCluster2, commonMuonPfoList);

    if (commonMuonPfoList.empty())
        return STATUS_CODE_NOT_FOUND;

    CartesianPointVector projectedPositions;
    StatusCode status(this->GetProjectedPositions(pCluster1, pCluster2, projectedPositions));

    if (status != STATUS_CODE_SUCCESS)
        return status;

    // Find all matched clusters (including unavailable)
    ClusterList matchedClusterList;
    this->CollectThirdViewClusters(pCluster1, pCluster2, projectedPositions, matchedClusterList);

    if (matchedClusterList.empty())
        return STATUS_CODE_NOT_FOUND;

    float reducedChiSquared(std::numeric_limits<float>::max());
    const Cluster *const pBestMatchedCluster =
        this->GetBestMatchedCluster(pCluster1, pCluster2, commonMuonPfoList, matchedClusterList, reducedChiSquared);

    //ATTN: Ignore if other clusters matches have more hits
    if (pBestMatchedCluster && (pBestMatchedCluster->IsAvailable()))
    {
        const unsigned int hitSum12(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits());
        const unsigned int hitSum13(pCluster1->GetNCaloHits() + pBestMatchedCluster->GetNCaloHits());

        if (hitSum13 > hitSum12)
            return STATUS_CODE_NOT_FOUND;

        const unsigned int hitSum23(pCluster2->GetNCaloHits() + pBestMatchedCluster->GetNCaloHits());

        if (hitSum23 > hitSum12)
            return STATUS_CODE_NOT_FOUND;
    }

    TwoViewXOverlap xOverlapObject(xMin1, xMax1, xMin2, xMax2);
    overlapResult = TwoViewDeltaRayOverlapResult(xOverlapObject, commonMuonPfoList, pBestMatchedCluster, matchedClusterList, reducedChiSquared);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::FindCommonMuonParents(const Cluster *const pCluster1, const Cluster *const pCluster2, PfoList &commonMuonPfoList) const
{
    ClusterList consideredClusters1, consideredClusters2;
    PfoList nearbyMuonPfos1, nearbyMuonPfos2;

    this->GetNearbyMuonPfos(pCluster1, consideredClusters1, nearbyMuonPfos1);

    if (nearbyMuonPfos1.empty())
        return;

    this->GetNearbyMuonPfos(pCluster2, consideredClusters2, nearbyMuonPfos2);

    if (nearbyMuonPfos2.empty())
        return;

    for (const ParticleFlowObject *const pNearbyMuon1 : nearbyMuonPfos1)
    {
        for (const ParticleFlowObject *const pNearbyMuon2 : nearbyMuonPfos2)
        {
            if (pNearbyMuon1 == pNearbyMuon2)
                commonMuonPfoList.push_back(pNearbyMuon1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::CollectThirdViewClusters(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const CartesianPointVector &projectedPositions, ClusterList &matchedClusters) const
{
    const ClusterList *pInputClusterList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));

    if (!pInputClusterList || pInputClusterList->empty())
        return;

    for (const Cluster *const pCluster : *pInputClusterList)
    {
        const float separation(LArMuonLeadingHelper::GetClosestDistance(pCluster, projectedPositions));

        if (separation > m_maxDistanceFromPrediction)
            continue;

        float reducedChiSquared(0.f);
        if (this->PerformThreeViewMatching(pCluster1, pCluster2, pCluster, reducedChiSquared) == STATUS_CODE_NOT_FOUND)
            continue;

        if (reducedChiSquared > m_maxGoodMatchReducedChiSquared)
            continue;

        matchedClusters.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TwoViewDeltaRayMatchingAlgorithm::GetBestMatchedCluster(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const PfoList &commonMuonPfoList, const ClusterList &matchedClusters, float &reducedChiSquared) const
{
    const Cluster *pBestMatchedCluster(nullptr);

    if (matchedClusters.empty())
        return pBestMatchedCluster;

    const HitType thirdViewHitType(LArClusterHelper::GetClusterHitType(matchedClusters.front()));
    ClusterList muonClusterList;

    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
        LArPfoHelper::GetClusters(pMuonPfo, thirdViewHitType, muonClusterList);

    unsigned int highestNHits(0);

    for (const Cluster *const pMatchedCluster : matchedClusters)
    {
        if (!pMatchedCluster->IsAvailable())
        {
            if (std::find(muonClusterList.begin(), muonClusterList.end(), pMatchedCluster) == muonClusterList.end())
                continue;
        }

        if (pMatchedCluster->GetNCaloHits() > highestNHits)
        {
            highestNHits = pMatchedCluster->GetNCaloHits();
            pBestMatchedCluster = pMatchedCluster;
        }
    }

    if (!pBestMatchedCluster)
        return pBestMatchedCluster;

    if (this->PerformThreeViewMatching(pCluster1, pCluster2, pBestMatchedCluster, reducedChiSquared) == STATUS_CODE_NOT_FOUND)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return pBestMatchedCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMatchingAlgorithm::CreatePfo(const MatrixType::Element &protoParticleElement)
{
    ProtoParticle protoParticle;

    protoParticle.m_clusterList.push_back(protoParticleElement.GetCluster1());
    protoParticle.m_clusterList.push_back(protoParticleElement.GetCluster2());

    const Cluster *const pBestMatchedCluster(protoParticleElement.GetOverlapResult().GetBestMatchedCluster());

    if (pBestMatchedCluster)
        this->FormThirdViewCluster(protoParticleElement, protoParticle);

    ProtoParticleVector protoParticleVector;
    protoParticleVector.emplace_back(protoParticle);

    return (this->CreatePfos(protoParticleVector));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::FormThirdViewCluster(const MatrixType::Element &element, ProtoParticle &protoParticle)
{
    const PfoList &commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());
    const Cluster *const pBestMatchedCluster(element.GetOverlapResult().GetBestMatchedCluster());
    const HitType thirdViewHitType(LArClusterHelper::GetClusterHitType(pBestMatchedCluster));
    const ParticleFlowObject *pMatchedMuonPfo(nullptr);

    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
    {
        ClusterList muonClusterList;
        LArPfoHelper::GetClusters(pMuonPfo, thirdViewHitType, muonClusterList);

        if (std::find(muonClusterList.begin(), muonClusterList.end(), pBestMatchedCluster) != muonClusterList.end())
            pMatchedMuonPfo = pMuonPfo;
    }

    const Cluster *pThirdViewCluster(pMatchedMuonPfo ? nullptr : pBestMatchedCluster);

    if (pMatchedMuonPfo)
    {
        CaloHitList deltaRayHitList;

        if (this->CollectHitsFromMuon(element.GetCluster1(), element.GetCluster2(), nullptr, pMatchedMuonPfo, m_minDistanceFromMuon,
                m_maxDistanceToCollected, deltaRayHitList) == STATUS_CODE_SUCCESS)
        {
            this->SplitMuonCluster(this->GetThirdViewClusterListName(), pBestMatchedCluster, deltaRayHitList, pThirdViewCluster);
            this->UpdateForThirdViewClusterModification(pBestMatchedCluster, true);
        }
        else
        {
            pThirdViewCluster = element.GetOverlapResult().GetBestMatchedAvailableCluster();
        }
    }

    if (!pThirdViewCluster)
        return;

    this->UpdateForThirdViewClusterModification(pThirdViewCluster, false);
    this->MergeThirdView(element, pThirdViewCluster);

    protoParticle.m_clusterList.push_back(pThirdViewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::MergeThirdView(const MatrixType::Element &element, const Cluster *const pSeedCluster)
{
    CaloHitList caloHitList1, caloHitList2;
    element.GetCluster1()->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
    element.GetCluster2()->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);

    // ATTN: Need copy as original will change throughout function
    ClusterList matchedClusters(element.GetOverlapResult().GetMatchedClusterList());

    ClusterSet checkedClusters;

    if (std::find(matchedClusters.begin(), matchedClusters.end(), pSeedCluster) != matchedClusters.end())
        checkedClusters.insert(pSeedCluster);

    while (checkedClusters.size() != matchedClusters.size())
    {
        const Cluster *pClusterToDelete(nullptr);
        unsigned int highestHit(0);

        for (const Cluster *const pMatchedCluster : matchedClusters)
        {
            if (checkedClusters.count(pMatchedCluster))
                continue;

            if (pMatchedCluster->GetNCaloHits() > highestHit)
            {
                pClusterToDelete = pMatchedCluster;
                highestHit = pMatchedCluster->GetNCaloHits();
            }
        }

        if (!pClusterToDelete)
            return;

        checkedClusters.insert(pClusterToDelete);

        if (!pClusterToDelete->IsAvailable())
            continue;

        CaloHitList caloHitList3;
        pSeedCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);
        pClusterToDelete->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);

        float reducedChiSquared(std::numeric_limits<float>::max());
        const StatusCode status(this->PerformThreeViewMatching(caloHitList1, caloHitList2, caloHitList3, reducedChiSquared));

        if (status == STATUS_CODE_NOT_FOUND)
            continue;

        if (reducedChiSquared > m_maxGoodMatchReducedChiSquared)
            continue;

        this->UpdateForThirdViewClusterModification(pClusterToDelete, false);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, this->GetThirdViewClusterListName()));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pClusterToDelete));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::UpdateForThirdViewClusterModification(const Cluster *const pModifiedCluster, const bool isMuon)
{
    auto &theMatrix(this->GetMatchingControl().GetOverlapMatrix());

    for (auto [pCluster1, overlapList] : theMatrix)
    {
        for (auto [pCluster2, overlapResult] : overlapList)
        {
            ClusterList matchedClusters(overlapResult.GetMatchedClusterList());

            auto matchedClustersIter(std::find(matchedClusters.begin(), matchedClusters.end(), pModifiedCluster));

            if (matchedClustersIter == matchedClusters.end())
                continue;

            float tempReducedChiSquared(std::numeric_limits<float>::max());

            if (isMuon)
                this->PerformThreeViewMatching(pCluster1, pCluster2, pModifiedCluster, tempReducedChiSquared);

            if (tempReducedChiSquared > m_maxGoodMatchReducedChiSquared)
                matchedClusters.erase(matchedClustersIter);

            float reducedChiSquared(std::numeric_limits<float>::max());
            const Cluster *const pBestMatchedCluster =
                this->GetBestMatchedCluster(pCluster1, pCluster2, overlapResult.GetCommonMuonPfoList(), matchedClusters, reducedChiSquared);

            TwoViewDeltaRayOverlapResult newOverlapResult(
                overlapResult.GetXOverlap(), overlapResult.GetCommonMuonPfoList(), pBestMatchedCluster, matchedClusters, reducedChiSquared);
            theMatrix.ReplaceOverlapResult(pCluster1, pCluster2, newOverlapResult);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    // Apply tools sequentially restarting if a change is made and ending if the tools finish or the restart limit is reached
    unsigned int repeatCounter(0);

    for (auto toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end();)
    {
        DeltaRayMatrixTool *const pTool(*toolIter);
        const bool repeatTools(pTool->Run(this, this->GetMatchingControl().GetOverlapMatrix()));

        toolIter = repeatTools ? m_algorithmToolVector.begin() : toolIter + 1;
        repeatCounter = repeatTools ? repeatCounter + 1 : repeatCounter;

        if (repeatCounter > m_nMaxMatrixToolRepeats)
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusterRebuilding", m_reclusteringAlgorithmName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "DeltaRayTools", algorithmToolVector));

    for (auto algorithmTool : algorithmToolVector)
    {
        DeltaRayMatrixTool *const pDeltaRayMatrixTool(dynamic_cast<DeltaRayMatrixTool *>(algorithmTool));

        if (!pDeltaRayMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayMatrixTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDistanceFromPrediction", m_maxDistanceFromPrediction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxGoodMatchReducedChiSquared", m_maxGoodMatchReducedChiSquared));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinDistanceFromMuon", m_minDistanceFromMuon));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToCollected", m_maxDistanceToCollected));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
