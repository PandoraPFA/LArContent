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
    m_minViewXSpan(0.f),
    m_minSeparation(1.f),
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

	usedKeyClusters.insert(pKeyCluster);

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

	this->SearchForMuonContamination(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::PassElementChecks(TensorType::Element &element, const HitType &hitType)
{
    const float viewXSpan(element.GetOverlapResult().GetViewXSpan(hitType));

    if (viewXSpan < m_minViewXSpan)
        return false;

    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
        return false;

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return false;

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(muonClusterList.front(), 10000, slidingFitPitch);

    const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList.front()));

    if (separation > m_minSeparation)
        return false;

    CartesianVector clusterAverageDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);

    const CartesianVector xAxis(1.f,0.f,0.f);

    if (!(((clusterAverageDirection.GetOpeningAngle(xAxis) * 180.f / 3.14) < 20.f) || (((180 - clusterAverageDirection.GetOpeningAngle(xAxis)) * 180.f / 3.14) < 20.f)))
    {
      std::cout << "not transverse" << std::endl; 
      return false;
    }

    CartesianVector deltaRayPoint(0.f,0.f,0.f), muonPoint(0.f,0.f,0.f);
    LArClusterHelper::GetClosestPositions(element.GetCluster(hitType), muonClusterList.front(), deltaRayPoint, muonPoint);

    return (this->ShouldSplitDeltaRay(muonClusterList.front(), element.GetCluster(hitType), muonPoint, slidingFitResult));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::SearchForMuonContamination(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
    ClusterSet modifiedClusters;
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (TensorType::Element &element : elementList)
    {
        for (const HitType &hitType : hitTypeVector)
	{
	    if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) || (modifiedClusters.count(element.GetClusterW())))
	        break;

	    if(!this->PassElementChecks(element, hitType))
	        continue;

	    const PfoList &commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

	    ClusterList muonClusterList;
	    LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);

	    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
	    const TwoDSlidingFitResult slidingFitResult(muonClusterList.front(), 20, slidingFitPitch);

	    // Is another element better?
	    bool isBestElement(true);
	    float bestChiSquared(element.GetOverlapResult().GetReducedChi2());
	    unsigned int highestHitNumber(element.GetClusterU()->GetNCaloHits() + element.GetClusterV()->GetNCaloHits() + element.GetClusterW()->GetNCaloHits()); 
	    for (TensorType::Element &otherElement : elementList)
	    {
	        if (otherElement.GetCluster(hitType) != element.GetCluster(hitType))
		    continue;

		if ((otherElement.GetClusterU() == element.GetClusterU()) && (otherElement.GetClusterV() == element.GetClusterV()) && (otherElement.GetClusterW() == element.GetClusterW()))
		    continue;

		unsigned int hitSum(otherElement.GetClusterU()->GetNCaloHits() + otherElement.GetClusterV()->GetNCaloHits() + otherElement.GetClusterW()->GetNCaloHits());
		if (hitSum < highestHitNumber)
		   continue;

		if (hitSum == highestHitNumber)
		{
		    if (otherElement.GetOverlapResult().GetReducedChi2() < bestChiSquared)
		        continue;
		}

		isBestElement = false;
	    }

	    if (!isBestElement)
	        continue;

	    // Attempt to pull delta ray hits out of cluster
            CaloHitList collectedHits;
	    this->CreateSeed(element, hitType, collectedHits);

	    if (collectedHits.empty())
	        continue;

            for (const CaloHit *const pCaloHit : collectedHits)
            {
	        const CartesianVector &position(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "seed", VIOLET, 2);
            }
	    PandoraMonitoringApi::ViewEvent(this->GetPandora());

	    if(!this->GrowSeed(element, hitType, collectedHits))
	      continue;

	    if (collectedHits.size() == element.GetCluster(hitType)->GetNCaloHits())
	        continue;

            for (const CaloHit *const pCaloHit : collectedHits)
            {
	        const CartesianVector &position(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "collected", VIOLET, 2);
            }

	    modifiedClusters.insert(element.GetCluster(hitType));

            std::cout << "xSpan: " << element.GetOverlapResult().GetViewXSpan(hitType) << std::endl;
	    ClusterList theCluster({element.GetCluster(hitType)});

	    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CLUSTER", BLUE);
	    PandoraMonitoringApi::ViewEvent(this->GetPandora());

	    this->SplitCluster(pAlgorithm, element, hitType, collectedHits);
	    changesMade = true;
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::ShouldSplitDeltaRay(const Cluster *const pMuonCluster, const Cluster *const pDeltaRayCluster, const CartesianVector &muonPosition,
    const TwoDSlidingFitResult &slidingFitResult) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(muonPosition, rL, rT);
    CartesianVector direction(0.f,0.f,0.f);
    slidingFitResult.GetGlobalFitDirection(rL, direction);

    CartesianVector minusPosition(muonPosition - (direction * 5.f));    
    CartesianVector plusPosition(muonPosition + (direction * 5.f));

    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &plusPosition, "plus", BLACK, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &minusPosition, "minus", BLACK, 2);

    CaloHitList minusMuonHits, minusDeltaRayHits;    
    CaloHitList plusMuonHits, plusDeltaRayHits;

    this->FindExtrapolatedHits(pMuonCluster, muonPosition, minusPosition, minusMuonHits);
    this->FindExtrapolatedHits(pMuonCluster, muonPosition, plusPosition, plusMuonHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonPosition, minusPosition, minusDeltaRayHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonPosition, plusPosition, plusDeltaRayHits);
    

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
   

    std::cout << "minusMuonHits.size(): " << minusMuonHits.size() << std::endl;
    std::cout << "plusMuonHits.size(): " << plusMuonHits.size() << std::endl;
    std::cout << "minusDeltaRayHits.size(): " << minusDeltaRayHits.size() << std::endl;
    std::cout << "plusDeltaRayHits.size(): " << plusDeltaRayHits.size() << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    if (minusMuonHits.empty() && plusMuonHits.empty())
        return false;

    if (minusDeltaRayHits.empty() && plusDeltaRayHits.empty())
        return false;

    /*
    // change this to be like, is the delta ray hits on the rhs close to that of the muons on the lhs?
    if ((minusMuonHits.size() < 3) && (plusDeltaRayHits.size() < 3) && (minusDeltaRayHits.size() > 2) && (plusMuonHits.size() > 2))
        return true;

    if ((plusMuonHits.size() < 3) && (minusDeltaRayHits.size() < 3) && (plusDeltaRayHits.size() > 2) && (minusMuonHits.size() > 2))
        return true; 
    

    if ((minusMuonHits.size() < 3) && (minusDeltaRayHits.size() > 2))
        return true;

    if ((plusMuonHits.size() < 3) && (plusDeltaRayHits.size() > 2))
        return true;   
    */

    if ((minusDeltaRayHits.size() > 2))
        return true;

    if ((plusDeltaRayHits.size() > 2))
        return true;   


    return false;
    
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

void CosmicRayRemovalTool::SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
  //delete the DR cluster and add in the new one (need to also delete and in the CR muon)
  pAlgorithm->UpdateUponDeletion(element.GetCluster(badHitType));

  ClusterList muonCluster;
  LArPfoHelper::GetClusters(element.GetOverlapResult().GetCommonMuonPfoList().front(), badHitType, muonCluster);
  pAlgorithm->UpdateUponDeletion(muonCluster.front());

  ClusterList oldDR({element.GetCluster(badHitType)});
  PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &oldDR, "old delta ray", BLUE);

  CaloHitList longClusterCaloHitList;
  element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(longClusterCaloHitList);
  for (const CaloHit *const pCaloHit : longClusterCaloHitList)
  {
    if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) == collectedHits.end())
    {
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, element.GetCluster(badHitType), pCaloHit));
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*pAlgorithm, muonCluster.front(), pCaloHit));
    }
  }

  ClusterList newDR({element.GetCluster(badHitType)});
  PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &newDR, "old delta ray", RED);
  PandoraMonitoringApi::ViewEvent(this->GetPandora());

  ClusterList newMuonCluster;
  LArPfoHelper::GetClusters(element.GetOverlapResult().GetCommonMuonPfoList().front(), badHitType, newMuonCluster);   

  ClusterVector clusterVector;
  clusterVector.push_back(newMuonCluster.front()); clusterVector.push_back(element.GetCluster(badHitType));//clusterVector.push_back(pDeltaRay);
  PfoVector pfoVector;
  pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front()); pfoVector.push_back(nullptr);

  pAlgorithm->UpdateForNewCluster(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::GrowSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), badHitType, muonClusterList);

    ClusterList muonClusterList1, muonClusterList2;

    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == badHitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == badHitType) || (hitType1 == hitType2))
                continue;

            LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType1, muonClusterList1);
            LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType2, muonClusterList2);

            if ((muonClusterList1.size() != 1) || (muonClusterList2.size() != 1))
            {
	      std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
	      return false;
	    }

	    break;
        }
	break;
    }
            
    CartesianPointVector muonProjectedPositions;
    this->ProjectPositions(muonClusterList1.front(), muonClusterList2.front(), muonProjectedPositions);

    for (const CartesianVector &position : muonProjectedPositions)
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", BLACK, 2);

    float projectedHitsFraction(static_cast<float>(muonProjectedPositions.size())  / muonClusterList.front()->GetNCaloHits());
    std::cout << "fraction of projected hits: " << static_cast<float>(muonProjectedPositions.size()) / muonClusterList.front()->GetNCaloHits() << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    CaloHitList longClusterCaloHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(longClusterCaloHitList);


    if (projectedHitsFraction < 0.8f)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
	const TwoDSlidingFitResult slidingFitResult(muonClusterList.front(), 40, slidingFitPitch);

	CartesianVector deltaRayPoint(0.f,0.f,0.f), muonPoint(0.f,0.f,0.f);
	LArClusterHelper::GetClosestPositions(element.GetCluster(badHitType), muonClusterList.front(), deltaRayPoint, muonPoint);

	CartesianVector closestProjectedPoint(this->GetClosestPosition(muonPoint, muonProjectedPositions, muonClusterList.front()));

	// want to return false so we stop the hit collection and know that it was incomplete
	if (closestProjectedPoint.GetMagnitude() < std::numeric_limits<float>::epsilon())
	   return false;

	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestProjectedPoint, "closest projection", RED, 2);
	PandoraMonitoringApi::ViewEvent(this->GetPandora());

	float rL(0.f), rT(0.f);
	CartesianVector clusterAverageDirection(0.f, 0.f, 0.f);
	slidingFitResult.GetLocalPosition(closestProjectedPoint, rL, rT);
	slidingFitResult.GetGlobalFitDirection(rL, clusterAverageDirection);

	std::cout << "impact param of muon point: " << clusterAverageDirection.GetCrossProduct(muonPoint - closestProjectedPoint).GetMagnitude() << std::endl;


	float zSpan(0.f);
        try
        {
	  float xSpanMin(muonPoint.GetX() - 1.f), xSpanMax(muonPoint.GetX() + 1.f), zSpanMin(0.f), zSpanMax(0.f);
	  muonClusterList.front()->GetClusterSpanZ(xSpanMin, xSpanMax, zSpanMin, zSpanMax);
	  zSpan = zSpanMax - zSpanMin;
	}
	catch (StatusCodeException &)
	{
	}


	std::cout << "cluster span:" << zSpan << std::endl;

	bool hitsAdded(true);
	while (hitsAdded)
        {
	    hitsAdded = false;

	    for (const CaloHit *const pCaloHit : longClusterCaloHitList)
            {
	        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
		  continue;

		//const float distanceToDeltaRayHitsSquared(this->GetClosestDistance(pCaloHit, collectedHits));
		const float distanceToMuonHitsSquared(clusterAverageDirection.GetCrossProduct(pCaloHit->GetPositionVector() - closestProjectedPoint).GetMagnitude());
		//float distanceToMuonHitsSquared(this->GetClosestDistance(pCaloHit, muonProjectedPositions));

	        if (distanceToMuonHitsSquared > 1.f)
	        {
	            collectedHits.push_back(pCaloHit);
		    hitsAdded = true;
		}
	    }
	}
    }
    else
    {
	bool hitsAdded(true);
	while (hitsAdded)
        {
	    hitsAdded = false;

	    for (const CaloHit *const pCaloHit : longClusterCaloHitList)
            {
	        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
		  continue;

		//const float distanceToDeltaRayHitsSquared(this->GetClosestDistance(pCaloHit, collectedHits));
		float distanceToMuonHitsSquared(this->GetClosestDistance(pCaloHit, muonProjectedPositions));

	        if (distanceToMuonHitsSquared > 0.5f)
	        {
	            collectedHits.push_back(pCaloHit);
		    hitsAdded = true;
		}
	    }
	}
    }

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

void CosmicRayRemovalTool::CreateSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
  CartesianPointVector projectedPositions, muonProjectedPositions;
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == badHitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == badHitType) || (hitType1 == hitType2))
                continue;


	    ClusterList muonClusterList1, muonClusterList2;
            LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType1, muonClusterList1);
            LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType2, muonClusterList2);

            if ((muonClusterList1.size() != 1) || (muonClusterList2.size() != 1))
            {
	      std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
	      return;
	    }

	    this->ProjectPositions(element.GetCluster(hitType1), element.GetCluster(hitType2), projectedPositions);
	    this->ProjectPositions(muonClusterList1.front(), muonClusterList2.front(), muonProjectedPositions);

	    break;
        }
	break;
    }

    CaloHitList longClusterCaloHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(longClusterCaloHitList);
    
    for (const CaloHit *const pCaloHit : longClusterCaloHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());

	float distanceToMuonHitsSquared(this->GetClosestDistance(pCaloHit, muonProjectedPositions));
	if (distanceToMuonHitsSquared < 1.f)
	  continue;

        for (const CartesianVector &projectedPosition : projectedPositions)
        {
            const float distanceSquared((position - projectedPosition).GetMagnitude());

            if (distanceSquared < 1.f)
	    {
                collectedHits.push_back(pCaloHit);
		break;
	    }
        }
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
  float shortestDistanceSquared(std::numeric_limits<float>::max());

  CartesianVector closestPoint(0.f,0.f,0.f);
  for (const CartesianVector &testPosition : cartesianPointVector)
  {
    if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > 0.5f)
      continue;

    float separationSquared((testPosition - referencePoint).GetMagnitude());

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

    
StatusCode CosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinViewXSpan", m_minViewXSpan));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapWindow", m_xOverlapWindow));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
