/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TwoViewCosmicRayRemovalTool.cc
 *
 *  @brief  Implementation of the two view cosmic removal tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewCosmicRayRemovalTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

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
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
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

        this->ExamineConnectedElements(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::ExamineConnectedElements(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList, bool &changesMade) const
{
    ClusterSet modifiedClusters, checkedClusters;
    
    for (const MatrixType::Element &element : elementList)
    {
        for (const HitType &hitType : pAlgorithm->GetHitTypeVector())
	    {
            const Cluster *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType));

            if (checkedClusters.count(pDeltaRayCluster))
                continue;
            
            if ((modifiedClusters.count(element.GetCluster1())) || (modifiedClusters.count(element.GetCluster2())))
                continue;

            if (!this->PassElementChecks(pAlgorithm, element, hitType))
                continue;

            if (!this->IsContaminated(pAlgorithm, element, hitType))
                continue;

            if (!this->IsBestElement(pAlgorithm, element, hitType, elementList))
                continue;

            checkedClusters.insert(pDeltaRayCluster);

	        // Attempt to pull delta ray hits out of cluster
            CaloHitList deltaRayHits;
            this->CreateSeed(pAlgorithm, element, hitType, deltaRayHits);

            if (deltaRayHits.empty())
                continue;

            // ATTN: If seed cannot be grown, abort
            CaloHitList deltaRayRemnantHits;
            if (this->GrowSeed(pAlgorithm, element, hitType, deltaRayHits, deltaRayRemnantHits) != STATUS_CODE_SUCCESS)
                continue;

            if (deltaRayHits.size() == pDeltaRayCluster->GetNCaloHits())
                continue;

            modifiedClusters.insert(pDeltaRayCluster);

            this->SplitDeltaRayCluster(pAlgorithm, element, hitType, deltaRayHits, deltaRayRemnantHits);
        }
    }

    changesMade = modifiedClusters.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::PassElementChecks(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType) const
{
    // ATTN: Avoid endpoints, topological michel reconstruction is very ambiguous
    if (this->IsMuonEndpoint(pAlgorithm, element, hitType))
        return false;

    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType));
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float separation(LArClusterHelper::GetClosestDistance(pDeltaRayCluster, pMuonCluster));
    
    if (separation > m_minSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsMuonEndpoint(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType) const
{
    for (const HitType &otherHitType : pAlgorithm->GetHitTypeVector())
    {
        if (otherHitType == hitType)
            continue;

        const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType));

        if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
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

bool TwoViewCosmicRayRemovalTool::IsContaminated(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType));
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
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

bool TwoViewCosmicRayRemovalTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, const CartesianVector &lineEnd, const float distanceToLine) const
{
    CartesianVector lineDirection(lineStart - lineEnd);
    lineDirection = lineDirection.GetUnitVector();
    
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
       return false;

    return true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
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

bool TwoViewCosmicRayRemovalTool::IsBestElement(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, 
    const HitType &hitType, const MatrixType::ElementList &elementList) const
{
    const  float chiSquared(element.GetOverlapResult().GetReducedChiSquared());
    const unsigned int hitSum(element.GetCluster1()->GetNCaloHits() + element.GetCluster2()->GetNCaloHits());
            
    for (const MatrixType::Element &testElement : elementList)
    {
        if (pAlgorithm->GetCluster(testElement, hitType) != pAlgorithm->GetCluster(element, hitType))
            continue;

        if ((testElement.GetCluster1() == element.GetCluster1()) && (testElement.GetCluster2() == element.GetCluster2()))
            continue;

        const  float testChiSquared(testElement.GetOverlapResult().GetReducedChiSquared());
        const unsigned int testHitSum(testElement.GetCluster1()->GetNCaloHits() + testElement.GetCluster2()->GetNCaloHits());

        if ((testHitSum < hitSum) || ((testHitSum == hitSum) && (testChiSquared > chiSquared)))
            continue;

        if (this->PassElementChecks(pAlgorithm, testElement, hitType))
            return false;
    }

    return true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::CreateSeed(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
    const HitType &hitType, CaloHitList &collectedHits) const
{
    // To avoid fluctuations, parameterise the muon track
    CartesianVector positionOnMuon(0.f, 0.f, 0.f), muonDirection(0.f, 0.f, 0.f);
    if (pAlgorithm->ParameteriseMuon(element.GetOverlapResult().GetCommonMuonPfoList().front(), pAlgorithm->GetCluster(element, hitType), 
        positionOnMuon, muonDirection) != STATUS_CODE_SUCCESS)
    {
        return;
    }

    CartesianPointVector deltaRayProjectedPositions;   
    if (this->ProjectDeltaRayPositions(pAlgorithm, element, hitType, deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CaloHitList deltaRayHitList;
    pAlgorithm->GetCluster(element, hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

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
    
StatusCode TwoViewCosmicRayRemovalTool::GrowSeed(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
    const HitType &hitType, CaloHitList &deltaRayHits, CaloHitList &remnantHits) const
{
    const Cluster *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType)), *pMuonCluster(nullptr);
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    // To avoid fluctuatuions, parameterise the muon track
    CartesianVector positionOnMuon(0.f, 0.f, 0.f), muonDirection(0.f, 0.f, 0.f);
    if (pAlgorithm->ParameteriseMuon(element.GetOverlapResult().GetCommonMuonPfoList().front(), pDeltaRayCluster, positionOnMuon, muonDirection) !=
        STATUS_CODE_SUCCESS)
    {
        return STATUS_CODE_NOT_FOUND;
    }

    // Identify delta ray hits
    this->CollectHitsFromDeltaRay(positionOnMuon, muonDirection, pDeltaRayCluster, m_maxDistanceToHit, true, deltaRayHits, deltaRayHits);

    // Identify remnant hits
    this->CollectHitsFromDeltaRay(positionOnMuon, muonDirection, pDeltaRayCluster, m_maxDistanceToHit, false, deltaRayHits, remnantHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::CollectHitsFromDeltaRay(const CartesianVector &positionOnMuon, const CartesianVector &muonDirection, const Cluster *const pDeltaRayCluster,
    const float &minDistanceFromMuon, const bool demandCloserToCollected, const CaloHitList &protectedHits, CaloHitList &collectedHits) const
{
    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);
    
    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

	    for (const CaloHit *const pCaloHit : deltaRayHitList)
        {
            if (std::find(protectedHits.begin(), protectedHits.end(), pCaloHit) != protectedHits.end())
                continue;            
            
            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const float distanceToCollectedHits(LArMuonLeadingHelper::GetClosestDistance(pCaloHit, collectedHits));
            const float distanceToMuonHits(muonDirection.GetCrossProduct(pCaloHit->GetPositionVector() - positionOnMuon).GetMagnitude());

            if ((distanceToMuonHits < minDistanceFromMuon) || (demandCloserToCollected && (distanceToCollectedHits > distanceToMuonHits)))
                continue;
                
            collectedHits.push_back(pCaloHit);
            hitsAdded = true;
        }
    }
} 

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::SplitDeltaRayCluster(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType,
    CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType)), *pMuonCluster(nullptr);    
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return;

    pAlgorithm->UpdateUponDeletion(pMuonCluster); pAlgorithm->UpdateUponDeletion(pDeltaRayCluster);

    std::string originalListName, fragmentListName;
    ClusterList originalClusterList(1, pDeltaRayCluster);
    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, originalClusterList, originalListName, fragmentListName));

    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    const Cluster *pDeltaRay(nullptr), *pDeltaRayRemnant(nullptr);

    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const bool isDeltaRay(std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end());
        const bool isDeltaRayRemnant(std::find(deltaRayRemnantHits.begin(), deltaRayRemnantHits.end(), pCaloHit) != deltaRayRemnantHits.end());
        const Cluster *&pCluster(isDeltaRay ? pDeltaRay : isDeltaRayRemnant ? pDeltaRayRemnant : pMuonCluster);
      
        if (!pCluster)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, parameters, pCluster));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*pAlgorithm, pCluster, pCaloHit));
        }
    }
    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*pAlgorithm, fragmentListName, originalListName));

    ClusterVector clusterVector; PfoVector pfoVector;
    if (pDeltaRayRemnant)
        this->FragmentRemnant(pAlgorithm, hitType, pMuonCluster, pDeltaRayRemnant, clusterVector, pfoVector);

    clusterVector.push_back(pMuonCluster); pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front());
    clusterVector.push_back(pDeltaRay); pfoVector.push_back(nullptr);

    pAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::FragmentRemnant(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const HitType &hitType, const Cluster *const pMuonCluster,
    const Cluster *const pDeltaRayRemnant, ClusterVector &clusterVector, PfoVector &pfoVector) const
{
    std::string caloHitListName(hitType == TPC_VIEW_U ? "CaloHitListU" : hitType == TPC_VIEW_V ? "CaloHitListV" : "CaloHitListW");
        
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*pAlgorithm, caloHitListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*pAlgorithm, pDeltaRayRemnant));

    const ClusterList *pClusterList(nullptr);
    std::string newClusterListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*pAlgorithm, pAlgorithm->GetClusteringAlgName(),
        pClusterList, newClusterListName));

    const ClusterList remnantClusters(pClusterList->begin(), pClusterList->end());
        
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, newClusterListName, pAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(hitType)));
        
    for (const Cluster *const pRemnant : remnantClusters)
    {
        if (pRemnant->GetNCaloHits() < m_minRemnantClusterSize)
        {
            if (LArClusterHelper::GetClosestDistance(pRemnant, pMuonCluster) < m_maxDistanceToTrack)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pMuonCluster, pRemnant));
                continue;
            }
        }

        clusterVector.push_back(pRemnant); pfoVector.push_back(nullptr);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::ProjectDeltaRayPositions(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
    const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const Cluster *pCluster1(nullptr), *pCluster2(nullptr);
    
    for (const HitType &hitType1 : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            pCluster1 = pAlgorithm->GetCluster(element, hitType1);
            pCluster2 = pAlgorithm->GetCluster(element, hitType2);
            
            break;
        }

        break;
    }

    if (!pCluster1 || !pCluster2)
        return STATUS_CODE_NOT_FOUND;

    return pAlgorithm->GetProjectedPositions(pCluster1, pCluster2, projectedPositions);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));   
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceToLine", m_distanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinContaminationLength", m_minContaminationLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceToHit", m_maxDistanceToHit));       

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinRemnantClusterSize", m_minRemnantClusterSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceToTrack", m_maxDistanceToTrack));      

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

