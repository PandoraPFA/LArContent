/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TwoViewCosmicRayRemovalTool.cc
 *
 *  @brief  Implementation of the two view cosmic removal tool class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewCosmicRayRemovalTool.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

TwoViewCosmicRayRemovalTool::TwoViewCosmicRayRemovalTool() :
    m_minSeparation(2.f),
    m_slidingFitWindow(10000),
    m_distanceToLine(0.5f),
    m_minContaminationLength(3.f),
    m_maxDistanceToHit(1.f),
    m_minRemnantClusterSize(3),
    m_maxDistanceToTrack(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    m_pParentAlgorithm = pAlgorithm;

    if (PandoraContentApi::GetSettings(*m_pParentAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    ClusterVector sortedKeyClusters;
    overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

    ClusterSet usedKeyClusters;
    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (usedKeyClusters.count(pKeyCluster))
            continue;

        MatrixType::ElementList elementList;
        overlapMatrix.GetConnectedElements(pKeyCluster, true, elementList);

        for (const MatrixType::Element &element : elementList)
            usedKeyClusters.insert(element.GetCluster1());

        const bool changesMadeInIteration = this->RemoveCosmicRayHits(elementList);

        changesMade = (changesMade ? changesMade : changesMadeInIteration);
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::RemoveCosmicRayHits(const MatrixType::ElementList &elementList) const
{
    ClusterSet modifiedClusters, checkedClusters;

    for (const MatrixType::Element &element : elementList)
    {
        for (const HitType hitType : m_pParentAlgorithm->GetHitTypeVector())
        {
            const Cluster *const pDeltaRayCluster(m_pParentAlgorithm->GetCluster(element, hitType));

            if (checkedClusters.count(pDeltaRayCluster))
                continue;

            // ATTN: The underlying matrix will update during this loop, do not proceed if element has been modified
            if ((modifiedClusters.count(element.GetCluster1())) || (modifiedClusters.count(element.GetCluster2())))
                continue;

            if (!this->PassElementChecks(element, hitType))
                continue;

            if (!this->IsContaminated(element, hitType))
                continue;

            if (!this->IsBestElement(element, hitType, elementList))
                continue;

            checkedClusters.insert(pDeltaRayCluster);

            // Attempt to pull delta ray hits out of cluster
            CaloHitList deltaRayHits;
            this->CreateSeed(element, hitType, deltaRayHits);

            if (deltaRayHits.empty())
                continue;

            // ATTN: If seed cannot be grown, abort
            CaloHitList deltaRayRemnantHits;
            if (this->GrowSeed(element, hitType, deltaRayHits, deltaRayRemnantHits) != STATUS_CODE_SUCCESS)
                continue;

            if (deltaRayHits.size() == pDeltaRayCluster->GetNCaloHits())
                continue;

            modifiedClusters.insert(pDeltaRayCluster);

            this->SplitDeltaRayCluster(element, hitType, deltaRayHits, deltaRayRemnantHits);
        }
    }

    return !modifiedClusters.empty();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::PassElementChecks(const MatrixType::Element &element, const HitType hitType) const
{
    // ATTN: Avoid endpoints, topological michel reconstruction is very ambiguous
    if (this->IsMuonEndpoint(element, hitType))
        return false;

    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(m_pParentAlgorithm->GetCluster(element, hitType));

    if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float separation(LArClusterHelper::GetClosestDistance(pDeltaRayCluster, pMuonCluster));

    return separation <= m_minSeparation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsMuonEndpoint(const MatrixType::Element &element, const HitType hitType) const
{
    for (const HitType &otherHitType : m_pParentAlgorithm->GetHitTypeVector())
    {
        if (otherHitType == hitType)
            continue;

        const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(m_pParentAlgorithm->GetCluster(element, hitType));

        if (m_pParentAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
            return false;

        float xMinDR(+std::numeric_limits<float>::max()), xMaxDR(-std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanX(xMinDR, xMaxDR);

        float xMinCR(-std::numeric_limits<float>::max()), xMaxCR(+std::numeric_limits<float>::max());
        pMuonCluster->GetClusterSpanX(xMinCR, xMaxCR);

        if ((xMinDR < xMinCR) || (xMaxDR > xMaxCR))
            return false;

        float zMinDR(+std::numeric_limits<float>::max()), zMaxDR(-std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanZ(xMinDR, xMaxDR, zMinDR, zMaxDR);

        float zMinCR(-std::numeric_limits<float>::max()), zMaxCR(+std::numeric_limits<float>::max());
        pMuonCluster->GetClusterSpanZ(xMinCR, xMaxCR, zMinCR, zMaxCR);

        if ((zMinDR < zMinCR) || (zMaxDR > zMaxCR))
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsContaminated(const MatrixType::Element &element, const HitType hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(m_pParentAlgorithm->GetCluster(element, hitType));

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

bool TwoViewCosmicRayRemovalTool::IsCloseToLine(
    const CartesianVector &hitPosition, const CartesianVector &lineStart, const CartesianVector &lineEnd, const float distanceToLine) const
{
    CartesianVector lineDirection(lineStart - lineEnd);
    lineDirection = lineDirection.GetUnitVector();

    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());

    if (transverseDistanceFromLine > distanceToLine)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsInLineSegment(
    const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float segmentBoundaryGradient = (-1.f) * (upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ()) / segmentBoundaryGradient + upperBoundary.GetX());
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ()) / segmentBoundaryGradient + lowerBoundary.GetX());

    if (std::fabs(xPointOnUpperLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;

    if (std::fabs(xPointOnLowerLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;

    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsBestElement(const MatrixType::Element &element, const HitType hitType, const MatrixType::ElementList &elementList) const
{
    const float chiSquared(element.GetOverlapResult().GetReducedChiSquared());
    const unsigned int hitSum(element.GetCluster1()->GetNCaloHits() + element.GetCluster2()->GetNCaloHits());

    for (const MatrixType::Element &testElement : elementList)
    {
        if (m_pParentAlgorithm->GetCluster(testElement, hitType) != m_pParentAlgorithm->GetCluster(element, hitType))
            continue;

        if ((testElement.GetCluster1() == element.GetCluster1()) && (testElement.GetCluster2() == element.GetCluster2()))
            continue;

        const float testChiSquared(testElement.GetOverlapResult().GetReducedChiSquared());
        const unsigned int testHitSum(testElement.GetCluster1()->GetNCaloHits() + testElement.GetCluster2()->GetNCaloHits());

        if ((testHitSum < hitSum) || ((testHitSum == hitSum) && (testChiSquared > chiSquared)))
            continue;

        if (this->PassElementChecks(testElement, hitType))
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::CreateSeed(const MatrixType::Element &element, const HitType hitType, CaloHitList &collectedHits) const
{
    // To avoid fluctuations, parameterise the muon track
    CartesianVector positionOnMuon(0.f, 0.f, 0.f), muonDirection(0.f, 0.f, 0.f);
    if (m_pParentAlgorithm->ParameteriseMuon(element.GetOverlapResult().GetCommonMuonPfoList().front(),
            m_pParentAlgorithm->GetCluster(element, hitType), positionOnMuon, muonDirection) != STATUS_CODE_SUCCESS)
        return;

    CartesianPointVector deltaRayProjectedPositions;
    if (this->ProjectDeltaRayPositions(element, hitType, deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CaloHitList deltaRayHitList;
    m_pParentAlgorithm->GetCluster(element, hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

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

StatusCode TwoViewCosmicRayRemovalTool::GrowSeed(
    const MatrixType::Element &element, const HitType hitType, CaloHitList &deltaRayHits, CaloHitList &remnantHits) const
{
    const Cluster *const pDeltaRayCluster(m_pParentAlgorithm->GetCluster(element, hitType)), *pMuonCluster(nullptr);

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

void TwoViewCosmicRayRemovalTool::CollectHitsFromDeltaRay(const CartesianVector &positionOnMuon, const CartesianVector &muonDirection,
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

void TwoViewCosmicRayRemovalTool::SplitDeltaRayCluster(
    const MatrixType::Element &element, const HitType hitType, CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *const pDeltaRayCluster(m_pParentAlgorithm->GetCluster(element, hitType)), *pMuonCluster(nullptr);

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

void TwoViewCosmicRayRemovalTool::ReclusterRemnant(const HitType hitType, const Cluster *const pMuonCluster,
    const Cluster *const pDeltaRayRemnant, ClusterVector &clusterVector, PfoVector &pfoVector) const
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

StatusCode TwoViewCosmicRayRemovalTool::ProjectDeltaRayPositions(
    const MatrixType::Element &element, const HitType hitType, CartesianPointVector &projectedPositions) const
{
    HitTypeVector hitTypes({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    hitTypes.erase(std::find(hitTypes.begin(), hitTypes.end(), hitType));

    const Cluster *const pCluster1(m_pParentAlgorithm->GetCluster(element, hitTypes[0]));
    const Cluster *const pCluster2(m_pParentAlgorithm->GetCluster(element, hitTypes[1]));

    if (!pCluster1 || !pCluster2)
        return STATUS_CODE_NOT_FOUND;

    return m_pParentAlgorithm->GetProjectedPositions(pCluster1, pCluster2, projectedPositions);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DistanceToLine", m_distanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinContaminationLength", m_minContaminationLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToHit", m_maxDistanceToHit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRemnantClusterSize", m_minRemnantClusterSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistanceToTrack", m_maxDistanceToTrack));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
