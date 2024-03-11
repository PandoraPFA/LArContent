/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CosmicRayRemovalTool.cc
 *
 *  @brief  Implementation of the cosmic ray removal tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayRemovalTool.h"

using namespace pandora;

namespace lar_content
{

CosmicRayRemovalTool::CosmicRayRemovalTool() :
    m_slidingFitWindow(10000),
    m_minContaminationLength(3.f),
    m_maxDistanceToHit(1.f),
    m_minRemnantClusterSize(3),
    m_maxDistanceToTrack(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
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

        const bool changesMadeInIteration = this->RemoveCosmicRayHits(elementList);

        changesMade = (changesMade ? changesMade : changesMadeInIteration);
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::RemoveCosmicRayHits(const TensorType::ElementList &elementList) const
{
    ClusterSet modifiedClusters, checkedClusters;

    for (const TensorType::Element &element : elementList)
    {
        for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const Cluster *const pDeltaRayCluster(element.GetCluster(hitType));

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
            this->CreateSeed(element, hitType, deltaRayHits);

            if (deltaRayHits.empty())
                continue;

            // ATTN: Abort if seed cannot be grown
            CaloHitList remnantHits;
            if (this->GrowSeed(element, hitType, deltaRayHits, remnantHits) != STATUS_CODE_SUCCESS)
                continue;

            // ATTN: Abort if reformed DR cluster
            if (deltaRayHits.size() == pDeltaRayCluster->GetNCaloHits())
                continue;

            modifiedClusters.insert(pDeltaRayCluster);

            this->SplitDeltaRayCluster(element, hitType, deltaRayHits, remnantHits);
        }
    }

    return !modifiedClusters.empty();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::PassElementChecks(const TensorType::Element &element, const HitType hitType) const
{
    // ATTN: Avoid endpoints, topological michel reconstruction is very ambiguous
    if (this->IsMuonEndpoint(element, true, hitType))
        return false;

    return RemovalBaseTool::PassElementChecks(element, hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::IsContaminated(const TensorType::Element &element, const HitType hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));

    if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, m_slidingFitWindow, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);

    // Find furthest point in DR cluster that lies on the projected muon track
    float furthestSeparation(0.f);
    CartesianVector extendedPoint(0.f, 0.f, 0.f);

    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        const float separation((position - muonVertex).GetMagnitude());

        if (separation > furthestSeparation)
        {
            if (this->IsCloseToLine(position, muonVertex, muonVertex + muonDirection, m_distanceToLine))
            {
                furthestSeparation = separation;
                extendedPoint = position;
            }
        }
    }

    // Check if delta ray has significant track contamination
    if (furthestSeparation < m_minContaminationLength)
        return false;

    // ATTN: Avoid cases where opening angle is small - it is easy to make mistakes in these instances
    CaloHitList muonHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonHitList);

    for (const CaloHit *const pCaloHit : muonHitList)
    {
        if (this->IsInLineSegment(deltaRayVertex, extendedPoint, pCaloHit->GetPositionVector()))
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::CreateSeed(const TensorType::Element &element, const HitType hitType, CaloHitList &collectedHits) const
{
    // To avoid fluctuations, parameterise the muon track
    CartesianVector positionOnMuon(0.f, 0.f, 0.f), muonDirection(0.f, 0.f, 0.f);
    if (m_pParentAlgorithm->ParameteriseMuon(element.GetOverlapResult().GetCommonMuonPfoList().front(), element.GetCluster(hitType),
            positionOnMuon, muonDirection) != STATUS_CODE_SUCCESS)
        return;

    CartesianPointVector deltaRayProjectedPositions;
    if (this->ProjectDeltaRayPositions(element, hitType, deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CaloHitList deltaRayHitList;
    element.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    CaloHitVector deltaRayHitVector(deltaRayHitList.begin(), deltaRayHitList.end());
    std::sort(deltaRayHitVector.begin(), deltaRayHitVector.end(), LArClusterHelper::SortHitsByPosition);

    for (const CaloHit *const pCaloHit : deltaRayHitVector)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());

        for (const CartesianVector &projectedPosition : deltaRayProjectedPositions)
        {
            const float distanceToProjectionSquared((position - projectedPosition).GetMagnitudeSquared());

            if (distanceToProjectionSquared < m_maxDistanceToHit * m_maxDistanceToHit)
            {
                // To prevent gappy seeds
                if ((!collectedHits.empty()) && (LArClusterHelper::GetClosestDistance(position, collectedHits) > m_maxDistanceToHit))
                    continue;

                const float distanceToMuonHitsSquared(muonDirection.GetCrossProduct(position - positionOnMuon).GetMagnitudeSquared());

                if (distanceToMuonHitsSquared < m_maxDistanceToHit * m_maxDistanceToHit)
                    continue;

                collectedHits.push_back(pCaloHit);

                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayRemovalTool::GrowSeed(const TensorType::Element &element, const HitType hitType, CaloHitList &deltaRayHits, CaloHitList &remnantHits) const
{
    const Cluster *const pDeltaRayCluster(element.GetCluster(hitType)), *pMuonCluster(nullptr);

    if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    // To avoid fluctuatuions, parameterise the muon track
    CartesianVector positionOnMuon(0.f, 0.f, 0.f), muonDirection(0.f, 0.f, 0.f);
    if (m_pParentAlgorithm->ParameteriseMuon(element.GetOverlapResult().GetCommonMuonPfoList().front(), pDeltaRayCluster, positionOnMuon,
            muonDirection) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    // Identify delta ray hits
    this->CollectHitsFromDeltaRay(positionOnMuon, muonDirection, pDeltaRayCluster, m_maxDistanceToHit, true, deltaRayHits, deltaRayHits);

    // Identify remnant hits
    this->CollectHitsFromDeltaRay(positionOnMuon, muonDirection, pDeltaRayCluster, m_maxDistanceToHit, false, deltaRayHits, remnantHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::CollectHitsFromDeltaRay(const CartesianVector &positionOnMuon, const CartesianVector &muonDirection,
    const Cluster *const pDeltaRayCluster, const float &minDistanceFromMuon, const bool demandCloserToCollected,
    const CaloHitList &protectedHits, CaloHitList &collectedHits) const
{
    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    bool hitsAdded(false);

    do
    {
        hitsAdded = false;

        for (const CaloHit *const pCaloHit : deltaRayHitList)
        {
            if (std::find(protectedHits.begin(), protectedHits.end(), pCaloHit) != protectedHits.end())
                continue;

            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const float distanceToCollectedHits(collectedHits.empty()
                    ? std::numeric_limits<float>::max()
                    : LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), collectedHits));
            const float distanceToMuonHits(muonDirection.GetCrossProduct(pCaloHit->GetPositionVector() - positionOnMuon).GetMagnitude());

            if ((distanceToMuonHits < minDistanceFromMuon) || (demandCloserToCollected && (distanceToCollectedHits > distanceToMuonHits)))
                continue;

            collectedHits.push_back(pCaloHit);
            hitsAdded = true;
        }
    } while (hitsAdded);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::SplitDeltaRayCluster(
    const TensorType::Element &element, const HitType hitType, CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *const pDeltaRayCluster(element.GetCluster(hitType)), *pMuonCluster(nullptr);

    if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return;

    m_pParentAlgorithm->UpdateUponDeletion(pMuonCluster);
    m_pParentAlgorithm->UpdateUponDeletion(pDeltaRayCluster);

    std::string originalListName, fragmentListName;
    ClusterList originalClusterList(1, pDeltaRayCluster);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::ReplaceCurrentList<Cluster>(*m_pParentAlgorithm, m_pParentAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::InitializeFragmentation(*m_pParentAlgorithm, originalClusterList, originalListName, fragmentListName));

    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    const Cluster *pDeltaRay(nullptr), *pDeltaRayRemnant(nullptr);

    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const bool isDeltaRay(std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end());
        const bool isDeltaRayRemnant(
            !isDeltaRay && (std::find(deltaRayRemnantHits.begin(), deltaRayRemnantHits.end(), pCaloHit) != deltaRayRemnantHits.end()));
        const Cluster *&pCluster(isDeltaRay ? pDeltaRay : isDeltaRayRemnant ? pDeltaRayRemnant : pMuonCluster);

        if (!pCluster)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*m_pParentAlgorithm, parameters, pCluster));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*m_pParentAlgorithm, pCluster, pCaloHit));
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*m_pParentAlgorithm, fragmentListName, originalListName));

    ClusterVector clusterVector;
    PfoVector pfoVector;
    if (pDeltaRayRemnant)
        this->ReclusterRemnant(hitType, pMuonCluster, pDeltaRayRemnant, clusterVector, pfoVector);

    clusterVector.push_back(pMuonCluster);
    pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front());
    clusterVector.push_back(pDeltaRay);
    pfoVector.push_back(nullptr);

    m_pParentAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::ReclusterRemnant(const HitType hitType, const Cluster *const pMuonCluster, const Cluster *const pDeltaRayRemnant,
    ClusterVector &clusterVector, PfoVector &pfoVector) const
{
    std::string caloHitListName(hitType == TPC_VIEW_U ? "CaloHitListU" : hitType == TPC_VIEW_V ? "CaloHitListV" : "CaloHitListW");

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*m_pParentAlgorithm, caloHitListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::ReplaceCurrentList<Cluster>(*m_pParentAlgorithm, m_pParentAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*m_pParentAlgorithm, pDeltaRayRemnant));

    const ClusterList *pClusterList(nullptr);
    std::string newClusterListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::RunClusteringAlgorithm(*m_pParentAlgorithm, m_pParentAlgorithm->GetClusteringAlgName(), pClusterList, newClusterListName));

    const ClusterList remnantClusters(pClusterList->begin(), pClusterList->end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::SaveList<Cluster>(*m_pParentAlgorithm, newClusterListName, m_pParentAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::ReplaceCurrentList<Cluster>(*m_pParentAlgorithm, m_pParentAlgorithm->GetClusterListName(hitType)));

    for (const Cluster *const pRemnant : remnantClusters)
    {
        if (pRemnant->GetNCaloHits() < m_minRemnantClusterSize)
        {
            if (LArClusterHelper::GetClosestDistance(pRemnant, pMuonCluster) < m_maxDistanceToTrack)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*m_pParentAlgorithm, pMuonCluster, pRemnant));
                continue;
            }
        }

        clusterVector.push_back(pRemnant);
        pfoVector.push_back(nullptr);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinContaminationLength", m_minContaminationLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToHit", m_maxDistanceToHit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRemnantClusterSize", m_minRemnantClusterSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToTrack", m_maxDistanceToTrack));

    return RemovalBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
