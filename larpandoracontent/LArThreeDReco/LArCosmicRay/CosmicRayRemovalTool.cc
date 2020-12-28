/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CosmicRayRemovalTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayRemovalTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

using namespace pandora;

namespace lar_content
{

CosmicRayRemovalTool::CosmicRayRemovalTool() :
    m_minSeparation(2.f),
    m_xOverlapWindow(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    ClusterSet usedKeyClusters;
    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (usedKeyClusters.count(pKeyCluster))
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        for (const TensorType::Element &element : elementList)
            usedKeyClusters.insert(element.GetCluster(TPC_VIEW_U));

        this->RemoveMuonHits(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::RemoveMuonHits(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &changesMade) const
{
    ClusterSet modifiedClusters, checkedClusters;
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    for (const TensorType::Element &element : elementList)
    {
        if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) || (modifiedClusters.count(element.GetClusterW())))
            continue;
        
        for (const HitType &hitType : hitTypeVector)
	    {
            if (checkedClusters.count(element.GetCluster(hitType)))
                continue;
            
            if (!this->PassElementChecks(element, hitType))
                continue;
            
            if (!this->IsContaminated(element, hitType))
                continue;

            if (!this->IsBestElement(element, hitType, elementList))
                continue;

            checkedClusters.insert(element.GetCluster(hitType));

	        // Attempt to pull delta ray hits out of cluster
            CaloHitList deltaRayHits;
            this->CreateSeed(element, hitType, deltaRayHits);

            if (deltaRayHits.empty())
                continue;

            // ATTN: If seed cannot be grown, abort
            CaloHitList deltaRayRemnantHits;
            if (this->GrowSeed(element, hitType, deltaRayHits, deltaRayRemnantHits) != STATUS_CODE_SUCCESS)
                continue;

            if (deltaRayHits.size() == element.GetCluster(hitType)->GetNCaloHits())
                continue;

            modifiedClusters.insert(element.GetCluster(hitType));

            this->SplitCluster(pAlgorithm, element, hitType, deltaRayHits, deltaRayRemnantHits);
        }
    }

    changesMade = modifiedClusters.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::PassElementChecks(const TensorType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), pMuonCluster));
    
    if (separation > m_minSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::IsContaminated(const TensorType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    CartesianVector deltaRayVertex(0.f,0.f,0.f), muonVertex(0.f,0.f,0.f);
    LArClusterHelper::GetClosestPositions(element.GetCluster(hitType), pMuonCluster, deltaRayVertex, muonVertex);    

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 10000, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const bool isTransverse(((muonDirection.GetOpeningAngle(xAxis) * 180.f / 3.14) < 20.f) || (((180 - muonDirection.GetOpeningAngle(xAxis)) * 180.f / 3.14) < 20.f));
    
    // If not transverse, check for significant muon contamination
    // TO DO - PERMORMACE SET THIS (COULD DO IF AND ELSE?) OR JUST USE THIS?
    if (!isTransverse)
    {
        CaloHitList deltaRayHitList;
        element.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

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
        if (furthestSeparation < 5.f)
            return false;

        // Rule out cases where muon follows DR
        CaloHitList muonHitList;
        pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonHitList);
        
        for (const CaloHit *const pCaloHit : muonHitList)
        {
            if (this->IsInLineSegment(deltaRayVertex, extendedPoint, pCaloHit->GetPositionVector()))
                return false;
        }
    }

    const CartesianVector minusPosition(muonVertex - (muonDirection * 5.f)), plusPosition(muonVertex + (muonDirection * 5.f));


    CaloHitList minusMuonHits, minusDeltaRayHits, plusMuonHits, plusDeltaRayHits;
    const Cluster *const pDeltaRayCluster(element.GetCluster(hitType));
    this->FindExtrapolatedHits(pMuonCluster, muonVertex, minusPosition, minusMuonHits);
    this->FindExtrapolatedHits(pMuonCluster, muonVertex, plusPosition, plusMuonHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, minusPosition, minusDeltaRayHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, plusPosition, plusDeltaRayHits);
   
    if (minusMuonHits.empty() && plusMuonHits.empty())
        return false;

    if (minusDeltaRayHits.empty() && plusDeltaRayHits.empty())
        return false;

    // OR does have DR majority?? i.e. is fraction over 1?
    if ((minusDeltaRayHits.size() > 2))
        return true;

    if ((plusDeltaRayHits.size() > 2))
        return true;   

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::IsBestElement(const TensorType::Element &element, const HitType &hitType, const TensorType::ElementList &elementList) const
{
    float chiSquared(element.GetOverlapResult().GetReducedChi2());
    unsigned int hitNumber(element.GetClusterU()->GetNCaloHits() + element.GetClusterV()->GetNCaloHits() + element.GetClusterW()->GetNCaloHits());
            
    for (const TensorType::Element &testElement : elementList)
    {
        if (testElement.GetCluster(hitType) != element.GetCluster(hitType))
            continue;

        if ((testElement.GetClusterU() == element.GetClusterU()) && (testElement.GetClusterV() == element.GetClusterV()) && (testElement.GetClusterW() == element.GetClusterW()))
            continue;

        const unsigned int hitSum(testElement.GetClusterU()->GetNCaloHits() + testElement.GetClusterV()->GetNCaloHits() + testElement.GetClusterW()->GetNCaloHits());

        if ((hitSum == hitNumber) && (testElement.GetOverlapResult().GetReducedChi2() < chiSquared))
            return false;
        
        if (hitSum > hitNumber)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::CreateSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
    CartesianPointVector muonProjectedPositions, deltaRayProjectedPositions;

    this->ProjectMuonPositions(element, badHitType, muonProjectedPositions);    
    this->ProjectDeltaRayPositions(element, badHitType, deltaRayProjectedPositions);

    CaloHitList deltaRayHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);
    
    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());

        for (const CartesianVector &projectedPosition : deltaRayProjectedPositions)
        {
            const float distanceToProjectionSquared((position - projectedPosition).GetMagnitudeSquared());

            if (distanceToProjectionSquared < 1.f * 1.f)
	        {
                // TO DO MAKE THIS SQUARED
                const float distanceToMuonHits(this->GetClosestDistance(pCaloHit, muonProjectedPositions));
        
                if (distanceToMuonHits < 1.f)
                    continue;
        
                collectedHits.push_back(pCaloHit);
                
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayRemovalTool::GrowSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    CartesianPointVector muonProjectedPositions;
    this->ProjectMuonPositions(element, badHitType, muonProjectedPositions);
    
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, badHitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    const float projectedHitsFraction(static_cast<float>(muonProjectedPositions.size()) / pMuonCluster->GetNCaloHits());

    CartesianVector muonDirection(0.f, 0.f, 0.f), positionOnMuon(0.f, 0.f, 0.f);
    if (projectedHitsFraction < 0.8f)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 40, slidingFitPitch);

        CartesianVector deltaRayVertex(0.f,0.f,0.f), muonVertex(0.f,0.f,0.f);
        LArClusterHelper::GetClosestPositions(element.GetCluster(badHitType), pMuonCluster, deltaRayVertex, muonVertex);

        positionOnMuon = this->GetClosestPosition(muonVertex, muonProjectedPositions, pMuonCluster);

	    if (positionOnMuon.GetMagnitude() < std::numeric_limits<float>::epsilon())
            return STATUS_CODE_NOT_FOUND;

        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(positionOnMuon, rL, rT);
        slidingFitResult.GetGlobalFitDirection(rL, muonDirection);
    }

    CaloHitList deltaRayHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

	    for (const CaloHit *const pCaloHit : deltaRayHitList)
        {
            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const float distanceToDeltaRayHits(this->GetClosestDistance(pCaloHit, collectedHits));
            const float distanceToMuonHits((projectedHitsFraction < 0.8f) ? muonDirection.GetCrossProduct(pCaloHit->GetPositionVector() - positionOnMuon).GetMagnitude() :
                this->GetClosestDistance(pCaloHit, muonProjectedPositions));

            if ((distanceToMuonHits > 1.f) && (distanceToDeltaRayHits < distanceToMuonHits))
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
                this->GetClosestDistance(pCaloHit, muonProjectedPositions));

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

void CosmicRayRemovalTool::SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &badHitType,
    CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, badHitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return;

    pAlgorithm->UpdateUponDeletion(pMuonCluster);    
    pAlgorithm->UpdateUponDeletion(element.GetCluster(badHitType));

    std::string originalListName, fragmentListName;
    ClusterList originalClusterList(1, element.GetCluster(badHitType));
    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(badHitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, originalClusterList, originalListName, fragmentListName));

    CaloHitList deltaRayHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

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
        this->FragmentRemnant(pAlgorithm, badHitType, pMuonCluster, pDeltaRayRemnant, clusterVector, pfoVector);
    
    clusterVector.push_back(pMuonCluster); pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front());
    clusterVector.push_back(pDeltaRay); pfoVector.push_back(nullptr);
    
    pAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::FragmentRemnant(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const HitType &badHitType, const Cluster *const pMuonCluster,
    const Cluster *const pDeltaRayRemnant, ClusterVector &clusterVector, PfoVector &pfoVector) const
{
    std::string caloHitListName(badHitType == TPC_VIEW_U ? "CaloHitListU" : badHitType == TPC_VIEW_V ? "CaloHitListV" : "CaloHitListW");
        
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*pAlgorithm, caloHitListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(badHitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*pAlgorithm, pDeltaRayRemnant));

    const ClusterList *pClusterList(nullptr);
    std::string newClusterListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*pAlgorithm, pAlgorithm->GetClusteringAlgName(),
        pClusterList, newClusterListName));

    const ClusterList remnantClusters(pClusterList->begin(), pClusterList->end());
        
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, newClusterListName, pAlgorithm->GetClusterListName(badHitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(badHitType)));
        
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
            
void CosmicRayRemovalTool::FindExtrapolatedHits(const Cluster *const pCluster, const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary,
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

bool CosmicRayRemovalTool::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
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

bool CosmicRayRemovalTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, const CartesianVector &lineEnd, const float distanceToLine) const
{
    CartesianVector lineDirection(lineStart - lineEnd);
    lineDirection = lineDirection.GetUnitVector();
    
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
       return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CosmicRayRemovalTool::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
{
  float shortestDistanceSquared(std::numeric_limits<float>::max());
  const CartesianVector referencePoint(pCaloHit->GetPositionVector());

  for (const CaloHit *const pTestCaloHit : caloHitList)
  {
    const CartesianVector &position(pTestCaloHit->GetPositionVector());
    float separationSquared((position - referencePoint).GetMagnitude());

    if (separationSquared < shortestDistanceSquared)
      shortestDistanceSquared = separationSquared;
  }

  return shortestDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CosmicRayRemovalTool::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
{
  float shortestDistanceSquared(std::numeric_limits<float>::max());
  const CartesianVector referencePoint(pCaloHit->GetPositionVector());

  for (const CartesianVector &testPosition : cartesianPointVector)
  {
    float separationSquared((testPosition - referencePoint).GetMagnitude());

    if (separationSquared < shortestDistanceSquared)
      shortestDistanceSquared = separationSquared;
  }

  return shortestDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::ProjectMuonPositions(const TensorType::Element &element, const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            const Cluster *pMuonCluster1(nullptr), *pMuonCluster2(nullptr);
            
            if (this->GetMuonCluster(element, hitType1, pMuonCluster1) != STATUS_CODE_SUCCESS)
                return;
                        
            if (this->GetMuonCluster(element, hitType2, pMuonCluster2) != STATUS_CODE_SUCCESS)
                return;

            this->ProjectPositions(pMuonCluster1, pMuonCluster2, projectedPositions);

            break;
        }
        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::ProjectDeltaRayPositions(const TensorType::Element &element, const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            this->ProjectPositions(element.GetCluster(hitType1), element.GetCluster(hitType2), projectedPositions);

            break;
        }
        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::ProjectPositions(const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianPointVector &projectedPositions) const
{    
    float xMin1(-std::numeric_limits<float>::max()), xMax1(+std::numeric_limits<float>::max());
    float xMin2(-std::numeric_limits<float>::max()), xMax2(+std::numeric_limits<float>::max());

    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, xMin2) - xPitch);
    const float xMax(std::min(xMax1, xMax2) + xPitch);
    const float xOverlap(xMax - xMin);

    // this has already been done...
    if (xOverlap < std::numeric_limits<float>::epsilon())
      return;
    
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    if (hitType1 == hitType2)
        return;

    const unsigned int nPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    // Projection into third view
    for (unsigned int n = 0; n < nPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMin1(0.f), zMin2(0.f), zMax1(0.f), zMax2(0.f);
            pCluster1->GetClusterSpanZ(xmin, xmax, zMin1, zMax1);
            pCluster2->GetClusterSpanZ(xmin, xmax, zMin2, zMax2);

            const float z1(0.5f * (zMin1 + zMax1));
            const float z2(0.5f * (zMin2 + zMax2));

            // could make use of the chi-squared?
            float chi2;
            CartesianVector projection(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, CartesianVector(x, 0.f, z1), CartesianVector(x, 0.f, z2), projection, chi2);

            projectedPositions.push_back(projection);
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (statusCodeException.GetStatusCode() != STATUS_CODE_NOT_FOUND)
                throw statusCodeException.GetStatusCode();

            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector CosmicRayRemovalTool::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector, const Cluster *const pCluster) const
{
    CartesianVector closestPoint(0.f,0.f,0.f);
    float shortestDistanceSquared(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > 0.5f)
            continue;

        const float separationSquared((testPosition - referencePoint).GetMagnitude());

        if (separationSquared > 5.f)
            continue;

        if (separationSquared < shortestDistanceSquared)
        {
            shortestDistanceSquared = separationSquared;
            closestPoint = testPosition;
        }
    }

    return closestPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayRemovalTool::GetMuonCluster(const TensorType::Element &element, const HitType &hitType, const Cluster *&pMuonCluster) const
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

StatusCode CosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapWindow", m_xOverlapWindow));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content


  //Pass Checks
  /*
  // trying to check for two connection points of the DR to the CR muon
  //const float sep(LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), muonClusterList.front()));
  if (sep < 2.f)
  {
    furthestSeparation = (pCaloHit->GetPositionVector() - muonPoint).GetMagnitude();
    correspondingPoint = pCaloHit->GetPositionVector();
  }
  */


    /*
    ClusterList jam({pDeltaRayCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &jam, "delta ray", BLACK);

    for (const CaloHit *const pCaloHit : minusMuonHits)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "muon", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : plusMuonHits)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "muon", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : minusDeltaRayHits)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "Dr", RED, 2);
    }

    for (const CaloHit *const pCaloHit : plusDeltaRayHits)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "Dr", RED, 2);
    }
   

    //std::cout << "minusMuonHits.size(): " << minusMuonHits.size() << std::endl;
    //std::cout << "plusMuonHits.size(): " << plusMuonHits.size() << std::endl;
    //std::cout << "minusDeltaRayHits.size(): " << minusDeltaRayHits.size() << std::endl;
    //std::cout << "plusDeltaRayHits.size(): " << plusDeltaRayHits.size() << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    // change this to be like, is the delta ray hits on the rhs close to that of the muons on the lhs?
    if ((minusMuonHits.size() < 3) && (plusDeltaRayHits.size() < 3) && (minusDeltaRayHits.size() > 2) && (plusMuonHits.size() > 2) && (minusDRTrackCount > 2))
        return true;

    if ((plusMuonHits.size() < 3) && (minusDeltaRayHits.size() < 3) && (plusDeltaRayHits.size() > 2) && (minusMuonHits.size() > 2) && (plusDRTrackCount > 2))
        return true; 
    

    if ((minusMuonHits.size() < 3) && (minusDeltaRayHits.size() > 2))
        return true;

    if ((plusMuonHits.size() < 3) && (plusDeltaRayHits.size() > 2))
        return true;   
    */
