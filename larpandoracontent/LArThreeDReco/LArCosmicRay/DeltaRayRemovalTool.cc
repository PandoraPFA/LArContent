/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.cc
 *
 *  @brief  Implementation of the delta ray removal tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.h"

using namespace pandora;

namespace lar_content
{

DeltaRayRemovalTool::DeltaRayRemovalTool() :
    m_slidingFitWindow(10000),
    m_minDeviationFromTransverse(0.35f),
    m_contaminationWindow(5.f),
    m_significantHitThreshold(3),
    m_minDistanceFromMuon(1.f),
    m_maxDistanceToCollected(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    m_pParentAlgorithm = pAlgorithm;

    if (PandoraContentApi::GetSettings(*m_pParentAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

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

        const bool changesMadeInIteration = this->RemoveDeltaRayHits(elementList);

        changesMade = (changesMade ? changesMade : changesMadeInIteration);
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::RemoveDeltaRayHits(const TensorType::ElementList &elementList) const
{
    ClusterSet modifiedClusters, checkedClusters;

    for (const TensorType::Element &element : elementList)
    {
        for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const Cluster *pDeltaRayCluster(element.GetCluster(hitType));
            const ParticleFlowObject *const pMuonPfo(element.GetOverlapResult().GetCommonMuonPfoList().front());

            if (checkedClusters.count(pDeltaRayCluster))
                continue;

            // ATTN: The underlying tensor will update during this loop, do not proceed if element has been modified
            if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) ||
                (modifiedClusters.count(element.GetClusterW())))
                continue;

            if (!this->PassElementChecks(element, hitType))
                continue;

            if (!this->IsContaminated(element, hitType))
                continue;

            if (!this->IsBestElement(element, hitType, elementList, modifiedClusters))
                continue;

            checkedClusters.insert(pDeltaRayCluster);

            CaloHitList deltaRayHits;
            if (m_pParentAlgorithm->CollectHitsFromMuon(nullptr, nullptr, pDeltaRayCluster, pMuonPfo, m_minDistanceFromMuon,
                    m_maxDistanceToCollected, deltaRayHits) != STATUS_CODE_SUCCESS)
                continue;

            modifiedClusters.insert(pDeltaRayCluster);

            this->SplitMuonCluster(element, hitType, deltaRayHits);
        }
    }

    return !modifiedClusters.empty();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::PassElementChecks(const TensorType::Element &element, const HitType hitType) const
{
    // ATTN: Avoid endpoints, topological michel reconstruction is very ambiguous
    if (this->IsMuonEndpoint(element, false))
        return false;

    return RemovalBaseTool::PassElementChecks(element, hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::IsContaminated(const TensorType::Element &element, const HitType hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));

    if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, m_slidingFitWindow, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const bool isTransverse((muonDirection.GetOpeningAngle(xAxis) < m_minDeviationFromTransverse) ||
        ((M_PI - muonDirection.GetOpeningAngle(xAxis)) < m_minDeviationFromTransverse));

    if (isTransverse)
        return false;

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);

    CaloHitList minusMuonHits, minusDeltaRayHits, plusMuonHits, plusDeltaRayHits;
    const CartesianVector minusPosition(muonVertex - (muonDirection * m_contaminationWindow));
    const CartesianVector plusPosition(muonVertex + (muonDirection * m_contaminationWindow));

    this->FindExtrapolatedHits(pMuonCluster, muonVertex, minusPosition, minusMuonHits);
    this->FindExtrapolatedHits(pMuonCluster, muonVertex, plusPosition, plusMuonHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, minusPosition, minusDeltaRayHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, plusPosition, plusDeltaRayHits);

    if ((minusMuonHits.size() < m_significantHitThreshold) && (plusMuonHits.size() < m_significantHitThreshold))
        return true;

    if ((minusMuonHits.size() < m_significantHitThreshold) && (minusDeltaRayHits.size() >= m_significantHitThreshold) &&
        (plusDeltaRayHits.size() < m_significantHitThreshold) && (plusMuonHits.size() >= m_significantHitThreshold))
    {
        return true;
    }

    if ((plusMuonHits.size() < m_significantHitThreshold) && (plusDeltaRayHits.size() >= m_significantHitThreshold) &&
        (minusDeltaRayHits.size() < m_significantHitThreshold) && (minusMuonHits.size() >= m_significantHitThreshold))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::SplitMuonCluster(const TensorType::Element &element, const HitType hitType, const CaloHitList &deltaRayHits) const
{
    const Cluster *pDeltaRayCluster(element.GetCluster(hitType)), *pMuonCluster(nullptr);

    if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_pParentAlgorithm->UpdateUponDeletion(pMuonCluster);
    m_pParentAlgorithm->UpdateUponDeletion(pDeltaRayCluster);

    m_pParentAlgorithm->SplitMuonCluster(m_pParentAlgorithm->GetClusterListName(hitType), pMuonCluster, deltaRayHits, pDeltaRayCluster);

    ClusterVector clusterVector;
    PfoVector pfoVector;
    clusterVector.push_back(pMuonCluster);
    pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front());
    clusterVector.push_back(pDeltaRayCluster);
    pfoVector.push_back(nullptr);
    m_pParentAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinDeviationFromTransverse", m_minDeviationFromTransverse));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ContaminationWindow", m_contaminationWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SignificantHitThreshold", m_significantHitThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinDistanceFromMuon", m_minDistanceFromMuon));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToCollected", m_maxDistanceToCollected));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
