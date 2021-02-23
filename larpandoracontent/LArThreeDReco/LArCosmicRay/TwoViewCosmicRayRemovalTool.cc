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
    m_minSeparation(2.f)
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

        this->RemoveMuonHits(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::RemoveMuonHits(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList, bool &changesMade) const
{
    const HitTypeVector hitTypeVector(pAlgorithm->GetHitTypeVector());
    
    ClusterSet modifiedClusters, checkedClusters;
    
    for (const MatrixType::Element &element : elementList)
    {
        for (const HitType &hitType : hitTypeVector)
	    {
	      if (checkedClusters.count(pAlgorithm->GetCluster(element, hitType)))
                continue;
            
            if ((modifiedClusters.count(element.GetCluster1())) || (modifiedClusters.count(element.GetCluster2())))
                continue;

            if (!this->PassElementChecks(pAlgorithm, element, hitType))
                continue;

            if (!this->IsContaminated(pAlgorithm, element, hitType))
                continue;

            if (!this->IsBestElement(pAlgorithm, element, hitType, elementList))
                continue;

            checkedClusters.insert(pAlgorithm->GetCluster(element, hitType));

	        // Attempt to pull delta ray hits out of cluster
            CaloHitList deltaRayHits;
            this->CreateSeed(pAlgorithm, element, hitType, deltaRayHits);

            if (deltaRayHits.empty())
                continue;

            // ATTN: If seed cannot be grown, abort
            CaloHitList deltaRayRemnantHits;
            if (this->GrowSeed(pAlgorithm, element, hitType, deltaRayHits, deltaRayRemnantHits) != STATUS_CODE_SUCCESS)
                continue;

            if (deltaRayHits.size() == pAlgorithm->GetCluster(element, hitType)->GetNCaloHits())
                continue;

            modifiedClusters.insert(pAlgorithm->GetCluster(element, hitType));

            this->SplitCluster(pAlgorithm, element, hitType, deltaRayHits, deltaRayRemnantHits);
        }
    }

    changesMade = modifiedClusters.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::PassElementChecks(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float separation(LArClusterHelper::GetClosestDistance(pAlgorithm->GetCluster(element, hitType), pMuonCluster));
    
    if (separation > m_minSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsContaminated(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType) const
{
    const HitTypeVector hitTypeVector(pAlgorithm->GetHitTypeVector());

    for (const HitType &otherHitType : hitTypeVector)
    {
        if (otherHitType == hitType)
            continue;

        const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType));
    
        if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
            return false;

        float xMinDR(+std::numeric_limits<float>::max()), xMaxDR(-std::numeric_limits<float>::max());
        float xMinCR(-std::numeric_limits<float>::max()), xMaxCR(+std::numeric_limits<float>::max());
        
        pDeltaRayCluster->GetClusterSpanX(xMinDR, xMaxDR);
        pMuonCluster->GetClusterSpanX(xMinCR, xMaxCR);

        if ((xMinDR < xMinCR) || (xMaxDR > xMaxCR))
            return false;
	
        float zMinDR(+std::numeric_limits<float>::max()), zMaxDR(-std::numeric_limits<float>::max());
        float zMinCR(-std::numeric_limits<float>::max()), zMaxCR(+std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanZ(xMinDR, xMaxDR, zMinDR, zMaxDR);
        pMuonCluster->GetClusterSpanZ(xMinCR, xMaxCR, zMinCR, zMaxCR);

        if ((zMinDR < zMinCR) || (zMaxDR > zMaxCR))
            return false;
    }

    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(pAlgorithm->GetCluster(element, hitType));
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);    

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 10000, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    float furthestSeparation(0.f); CartesianVector extendedPoint(0.f,0.f,0.f);
        
    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        const float separation((position - muonVertex).GetMagnitude());
            
        if (separation > furthestSeparation)
	    {
            if (this->IsCloseToLine(position, muonVertex, muonVertex + muonDirection, 0.5f))
            {
                furthestSeparation = separation;
                extendedPoint = position;
            }
        }
    }


    // Check if significant
    if (furthestSeparation < 3.f)
        return false;

    // Rule out cases where muon follows DR
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

void TwoViewCosmicRayRemovalTool::CreateSeed(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
    const HitType &hitType, CaloHitList &collectedHits) const
{
    CartesianPointVector muonProjectedPositions;
    if (pAlgorithm->ProjectMuonPositions(hitType, element.GetOverlapResult().GetCommonMuonPfoList().front(), muonProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CartesianPointVector deltaRayProjectedPositions;   
    if (this->ProjectDeltaRayPositions(pAlgorithm, element, hitType, deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CaloHitList deltaRayHitList;
    pAlgorithm->GetCluster(element, hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);
    
    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());

        for (const CartesianVector &projectedPosition : deltaRayProjectedPositions)
        {
            const float distanceToProjectionSquared((position - projectedPosition).GetMagnitudeSquared());

            if (distanceToProjectionSquared < 1.f * 1.f)
	        {
                // TO DO MAKE THIS SQUARED
                const float distanceToMuonHits(LArMuonLeadingHelper::GetClosestDistance(pCaloHit, muonProjectedPositions));
        
                if (distanceToMuonHits < 1.f)
                    continue;
        
                collectedHits.push_back(pCaloHit);
                
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::GrowSeed(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
    const HitType &hitType, CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    CartesianPointVector muonProjectedPositions;
    if (pAlgorithm->ProjectMuonPositions(hitType, element.GetOverlapResult().GetCommonMuonPfoList().front(), muonProjectedPositions) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    const float projectedHitsFraction(static_cast<float>(muonProjectedPositions.size()) / pMuonCluster->GetNCaloHits());

    CartesianVector muonDirection(0.f, 0.f, 0.f), positionOnMuon(0.f, 0.f, 0.f);
    if (projectedHitsFraction < 0.8f)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 40, slidingFitPitch);

        CartesianVector deltaRayVertex(0.f,0.f,0.f), muonVertex(0.f,0.f,0.f);
        LArClusterHelper::GetClosestPositions(pAlgorithm->GetCluster(element, hitType), pMuonCluster, deltaRayVertex, muonVertex);

        positionOnMuon = LArMuonLeadingHelper::GetClosestPosition(muonVertex, muonProjectedPositions, pMuonCluster);

	    if (positionOnMuon.GetMagnitude() < std::numeric_limits<float>::epsilon())
            return STATUS_CODE_NOT_FOUND;

        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(positionOnMuon, rL, rT);
        slidingFitResult.GetGlobalFitDirection(rL, muonDirection);
    }

    CaloHitList deltaRayHitList;
    pAlgorithm->GetCluster(element, hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

	    for (const CaloHit *const pCaloHit : deltaRayHitList)
        {
            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const float distanceToDeltaRayHits(LArMuonLeadingHelper::GetClosestDistance(pCaloHit, collectedHits));
            const float distanceToMuonHits((projectedHitsFraction < 0.8f) ? muonDirection.GetCrossProduct(pCaloHit->GetPositionVector() - positionOnMuon).GetMagnitude() :
                LArMuonLeadingHelper::GetClosestDistance(pCaloHit, muonProjectedPositions));

            if ((distanceToMuonHits > 0.5f) && (distanceToDeltaRayHits < distanceToMuonHits))
	        {
                collectedHits.push_back(pCaloHit);
                hitsAdded = true;
            }
        }
    }

    hitsAdded = true;
	while (hitsAdded)
    {
	    hitsAdded = false;

	    for (const CaloHit *const pCaloHit : deltaRayHitList)
        {
	        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

	        if (std::find(deltaRayRemnantHits.begin(), deltaRayRemnantHits.end(), pCaloHit) != deltaRayRemnantHits.end())
                continue;

            const float distanceToMuonHits((projectedHitsFraction < 0.8f) ? muonDirection.GetCrossProduct(pCaloHit->GetPositionVector() - positionOnMuon).GetMagnitude() :
                LArMuonLeadingHelper::GetClosestDistance(pCaloHit, muonProjectedPositions));

	        if (distanceToMuonHits > 1.f)
	        {
	            deltaRayRemnantHits.push_back(pCaloHit);
                hitsAdded = true;
            }
	    }
	}

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewCosmicRayRemovalTool::SplitCluster(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType,
    CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return;

    pAlgorithm->UpdateUponDeletion(pMuonCluster);    
    pAlgorithm->UpdateUponDeletion(pAlgorithm->GetCluster(element, hitType));

    std::string originalListName, fragmentListName;
    ClusterList originalClusterList(1, pAlgorithm->GetCluster(element, hitType));
    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, originalClusterList, originalListName, fragmentListName));

    CaloHitList deltaRayHitList;
    pAlgorithm->GetCluster(element, hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

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
        if (pRemnant->GetNCaloHits() < 3)
        {
            if (LArClusterHelper::GetClosestDistance(pRemnant, pMuonCluster) < 2.f)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pMuonCluster, pRemnant));
                continue;
            }
        }

        clusterVector.push_back(pRemnant); pfoVector.push_back(nullptr);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::GetMuonCluster(const MatrixType::Element &element, const HitType &hitType, const Cluster *&pMuonCluster) const
{
    const PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    pMuonCluster = muonClusterList.front();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewCosmicRayRemovalTool::IsBestElement(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const HitType &hitType, const MatrixType::ElementList &elementList) const
{
    float chiSquared(element.GetOverlapResult().GetReducedChiSquared());
    unsigned int hitNumber(element.GetCluster1()->GetNCaloHits() + element.GetCluster2()->GetNCaloHits()); //maybe add other in?
            
    for (const MatrixType::Element &testElement : elementList)
    {
        if (pAlgorithm->GetCluster(testElement, hitType) != pAlgorithm->GetCluster(element, hitType))
            continue;

        if ((testElement.GetCluster1() == element.GetCluster1()) && (testElement.GetCluster2() == element.GetCluster2()))
            continue;

        const unsigned int hitSum(testElement.GetCluster1()->GetNCaloHits() + testElement.GetCluster2()->GetNCaloHits());

        if ((hitSum == hitNumber) && (testElement.GetOverlapResult().GetReducedChiSquared() < chiSquared))
            return false;
        
        if (hitSum > hitNumber)
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

void TwoViewCosmicRayRemovalTool::FindExtrapolatedHits(const Cluster *const pCluster, const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary,
    CaloHitList &collectedHits) const
{
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (!this->IsInLineSegment(lowerBoundary, upperBoundary, pCaloHit->GetPositionVector()))
            continue;

        if (!this->IsCloseToLine(pCaloHit->GetPositionVector(), lowerBoundary, upperBoundary, 0.5))
            continue;

        collectedHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::ProjectDeltaRayPositions(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
    const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    const Cluster *pCluster1(nullptr), *pCluster2(nullptr);
    
    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            pCluster1 = pAlgorithm->GetCluster(element, hitType1);
            pCluster2 = pAlgorithm->GetCluster(element, hitType2);
            
            break;
        }
        break;
    }

    //ATTN: returned cluster may be a nullptr if no available best matched cluster
    if (!pCluster1 || !pCluster2)
        return STATUS_CODE_NOT_FOUND;

    return pAlgorithm->GetProjectedPositions(pCluster1, pCluster2, projectedPositions);
}



//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewCosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

