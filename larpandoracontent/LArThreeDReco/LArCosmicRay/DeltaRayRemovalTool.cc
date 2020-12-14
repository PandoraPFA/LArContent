/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/DeltaRayRemovalTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

using namespace pandora;

namespace lar_content
{

DeltaRayRemovalTool::DeltaRayRemovalTool() :
    m_minSeparation(1.f),
    m_xOverlapWindow(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
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

	this->SearchForDeltaRayContamination(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::PassElementChecks(TensorType::Element &element, const HitType &hitType)
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
    {
      std::cout << "common muon issue" << std::endl;
        return false;
    }

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return false;

    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    for (const HitType &jam : hitTypeVector)
    {
        ClusterList hitMuonClusterList;
	LArPfoHelper::GetClusters(commonMuonPfoList.front(), jam, hitMuonClusterList);
            
	if (hitMuonClusterList.size() != 1)
            return false;

        float xMinDR(+std::numeric_limits<float>::max()), xMaxDR(-std::numeric_limits<float>::max());
	float xMinCR(-std::numeric_limits<float>::max()), xMaxCR(+std::numeric_limits<float>::max());
	element.GetCluster(jam)->GetClusterSpanX(xMinDR, xMaxDR);
	hitMuonClusterList.front()->GetClusterSpanX(xMinCR, xMaxCR);

        if ((xMinDR < xMinCR) || (xMaxDR > xMaxCR))
        {
	  std::cout << "DR lies outside CR" << std::endl;
	  return false;
	}
	
        float zMinDR(+std::numeric_limits<float>::max()), zMaxDR(-std::numeric_limits<float>::max());
	float zMinCR(-std::numeric_limits<float>::max()), zMaxCR(+std::numeric_limits<float>::max());
	element.GetCluster(jam)->GetClusterSpanZ(xMinDR, xMaxDR, zMinDR, zMaxDR);
	hitMuonClusterList.front()->GetClusterSpanZ(xMinCR, xMaxCR, zMinCR, zMaxCR);

        if ((zMinDR < zMinCR) || (zMaxDR > zMaxCR))
        {
	  std::cout << "DR lies outside CR" << std::endl;
	  return false;
	}
    }

    const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList.front()));

    if (separation > m_minSeparation)
    {
      std::cout << "separation issue" << std::endl;
        return false;
    }

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(muonClusterList.front(), 100000, slidingFitPitch);

    CartesianVector clusterAverageDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);

    const CartesianVector xAxis(1.f,0.f,0.f);

    if (((clusterAverageDirection.GetOpeningAngle(xAxis) * 180.f / 3.14) < 20.f) || (((180 - clusterAverageDirection.GetOpeningAngle(xAxis)) * 180.f / 3.14) < 20.f))
    {
      std::cout << "too transverse" << std::endl; 
      return false;
    }

    CartesianVector deltaRayPoint(0.f,0.f,0.f), muonPoint(0.f,0.f,0.f);
    LArClusterHelper::GetClosestPositions(element.GetCluster(hitType), muonClusterList.front(), deltaRayPoint, muonPoint);

    return (this->ShouldSplitCosmicRay(muonClusterList.front(), element.GetCluster(hitType), muonPoint, slidingFitResult));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::SearchForDeltaRayContamination(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
    ClusterSet modifiedClusters;
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    std::cout << "trying to remove DR hits from CRs" << std::endl;

    for (TensorType::Element &element : elementList)
    {
        for (const HitType &hitType : hitTypeVector)
	{
	    if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) || (modifiedClusters.count(element.GetClusterW())))
	        break;

	    if (!this->PassElementChecks(element, hitType))
	        continue;

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
	    {
	        std::cout << "not the best elemet - search again" << std::endl;
	        continue;
	    }

	    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());
	    
	    ClusterList muonClusterList;
	    LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);

	    ClusterList theCluster({element.GetCluster(hitType)});
	    ClusterList muonCluster({muonClusterList.front()});
	    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CLUSTER", RED);
	    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &muonCluster, "CLUSTER", BLUE);
	    PandoraMonitoringApi::ViewEvent(this->GetPandora()); 
	    

	    // Attempt to pull delta ray hits out of muon cluster
            CaloHitList collectedHits;
	    this->CreateSeed(element, hitType, collectedHits);

	    if (collectedHits.empty())
	    {
	        std::cout << "collected hits are empty" << std::endl;
	        continue;
	    }

            for (const CaloHit *const pCaloHit : collectedHits)
            {
	        const CartesianVector &position(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "collected", VIOLET, 2);
            }
	    PandoraMonitoringApi::ViewEvent(this->GetPandora());

	    modifiedClusters.insert(element.GetCluster(hitType));

	    this->SplitCluster(pAlgorithm, element, hitType, collectedHits);
	    changesMade = true;
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::ShouldSplitCosmicRay(const Cluster *const pMuonCluster, const Cluster *const pDeltaRayCluster, const CartesianVector &muonPosition,
    const TwoDSlidingFitResult &slidingFitResult) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(muonPosition, rL, rT);
    CartesianVector direction(0.f,0.f,0.f);
    slidingFitResult.GetGlobalFitDirection(rL, direction);

    CartesianVector minusPosition(muonPosition - (direction * 5.f));    
    CartesianVector plusPosition(muonPosition + (direction * 5.f));

    CaloHitList minusMuonHits, minusDeltaRayHits;    
    CaloHitList plusMuonHits, plusDeltaRayHits;

    this->FindExtrapolatedHits(pMuonCluster, muonPosition, minusPosition, minusMuonHits);
    this->FindExtrapolatedHits(pMuonCluster, muonPosition, plusPosition, plusMuonHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonPosition, minusPosition, minusDeltaRayHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonPosition, plusPosition, plusDeltaRayHits);

    /*
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
    */
    std::cout << "minusMuonHits.size(): " << minusMuonHits.size() << std::endl;
    std::cout << "plusMuonHits.size(): " << plusMuonHits.size() << std::endl;
    std::cout << "minusDeltaRayHits.size(): " << minusDeltaRayHits.size() << std::endl;
    std::cout << "plusDeltaRayHits.size(): " << plusDeltaRayHits.size() << std::endl;

    // Are there no muon hits either side?
    if ((minusMuonHits.size() < 3) && (plusMuonHits.size() < 3))
        return true;

    // change this to be like, is the delta ray hits on the rhs close to that of the muons on the lhs?
    if ((minusMuonHits.size() < 3) && (plusDeltaRayHits.size() < 3) && (minusDeltaRayHits.size() > 2) && (plusMuonHits.size() > 2))
        return true;

    if ((plusMuonHits.size() < 3) && (minusDeltaRayHits.size() < 3) && (plusDeltaRayHits.size() > 2) && (minusMuonHits.size() > 2))
        return true;  



    return false;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::FindExtrapolatedHits(const Cluster *const pCluster, const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary,
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

bool DeltaRayRemovalTool::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
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

bool DeltaRayRemovalTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, const CartesianVector &lineEnd, const float distanceToLine) const
{
    CartesianVector lineDirection(lineStart - lineEnd);
    lineDirection = lineDirection.GetUnitVector();
    
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
       return false;

    return true;
}
//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::CreateSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    CartesianPointVector deltaRayProjectedPositions, muonProjectedPositions;

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

	    this->ProjectPositions(element.GetCluster(hitType1), element.GetCluster(hitType2), deltaRayProjectedPositions);
	    this->ProjectPositions(muonClusterList1.front(), muonClusterList2.front(), muonProjectedPositions);

	    break;
        }
	break;
    }

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), badHitType, muonClusterList);

    if (muonClusterList.size() != 1)
    {
        std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
	return;
    }

    for (const CartesianVector &position : deltaRayProjectedPositions)
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "DR", RED, 2);

    for (const CartesianVector &position : muonProjectedPositions)
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "CR", BLUE, 2);

    float projectedHitsFraction(static_cast<float>(muonProjectedPositions.size())  / muonClusterList.front()->GetNCaloHits());
    std::cout << "fraction of projected hits: " << static_cast<float>(muonProjectedPositions.size()) / muonClusterList.front()->GetNCaloHits() << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    CaloHitList muonCaloHitList;
    muonClusterList.front()->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);

    if (projectedHitsFraction < 1.f)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
	const TwoDSlidingFitResult slidingFitResult(muonClusterList.front(), 20, slidingFitPitch);

	CartesianVector clusterAverageDirection(0.f, 0.f, 0.f);
	//slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);

	CartesianVector deltaRayPoint(0.f,0.f,0.f), muonPoint(0.f,0.f,0.f);
	LArClusterHelper::GetClosestPositions(element.GetCluster(badHitType), muonClusterList.front(), deltaRayPoint, muonPoint);

	CartesianVector closestProjectedPoint(this->GetClosestPosition(muonPoint, muonProjectedPositions, muonClusterList.front()));

	if (closestProjectedPoint.GetMagnitude() < std::numeric_limits<float>::epsilon())
	   return;

	float rL(0.f), rT(0.f);
	slidingFitResult.GetLocalPosition(closestProjectedPoint, rL, rT);
	slidingFitResult.GetGlobalFitDirection(rL, clusterAverageDirection);

	bool hitsAdded(true);
	while (hitsAdded)
        {
	    hitsAdded = false;

            for (const CaloHit *const pCaloHit : muonCaloHitList)
            {
	        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
	            continue;

		//std::min(this->GetClosestDistance(pCaloHit, deltaRayProjectedPositions), 
		float distanceToDeltaRayHitsSquared(std::min(this->GetClosestDistance(pCaloHit, collectedHits), 
							     LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), element.GetCluster(badHitType))));
		const float distanceToMuonHitsSquared(clusterAverageDirection.GetCrossProduct(pCaloHit->GetPositionVector() - closestProjectedPoint).GetMagnitude());

	        if ((std::fabs(distanceToMuonHitsSquared - distanceToDeltaRayHitsSquared) > std::numeric_limits<float>::epsilon()) && (distanceToDeltaRayHitsSquared < distanceToMuonHitsSquared)
		    && (distanceToMuonHitsSquared > 1.f) && (distanceToDeltaRayHitsSquared < 1.f))
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

            for (const CaloHit *const pCaloHit : muonCaloHitList)
            {
	        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
	            continue;

		float distanceToDeltaRayHitsSquared(std::min(this->GetClosestDistance(pCaloHit, collectedHits), 
							     LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), element.GetCluster(badHitType))));
		float distanceToMuonHitsSquared(this->GetClosestDistance(pCaloHit, muonProjectedPositions));

		if ((std::fabs(distanceToMuonHitsSquared - distanceToDeltaRayHitsSquared) > std::numeric_limits<float>::epsilon()) && (distanceToDeltaRayHitsSquared < distanceToMuonHitsSquared)
		    && (distanceToMuonHitsSquared > 1.f) && (distanceToDeltaRayHitsSquared < 1.f))
	        {
	            collectedHits.push_back(pCaloHit);
		    hitsAdded = true;
		}
	    }
	}
    }

    // Catch if delta ray has travelled along muon
    std::cout << "collected hit fraction: " << (static_cast<float>(collectedHits.size()) / muonCaloHitList.size()) << std::endl;
    if ((static_cast<float>(collectedHits.size()) / muonCaloHitList.size()) > 0.2)
      collectedHits.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
  //delete the DR cluster and add in the new one (need to also delete and in the CR muon)
  pAlgorithm->UpdateUponDeletion(element.GetCluster(badHitType));

  ClusterList muonCluster;
  LArPfoHelper::GetClusters(element.GetOverlapResult().GetCommonMuonPfoList().front(), badHitType, muonCluster);
  pAlgorithm->UpdateUponDeletion(muonCluster.front());

  ClusterList oldDR({element.GetCluster(badHitType)});
  PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &oldDR, "old delta ray", BLUE);
  PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &muonCluster, "old MUON ray", BLACK);

  CaloHitList muonCaloHitList;
  muonCluster.front()->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);
  for (const CaloHit *const pCaloHit : muonCaloHitList)
  {
    if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
    {
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, muonCluster.front(), pCaloHit));
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*pAlgorithm, element.GetCluster(badHitType), pCaloHit));
    }
  }

  ClusterList newMuonCluster;
  LArPfoHelper::GetClusters(element.GetOverlapResult().GetCommonMuonPfoList().front(), badHitType, newMuonCluster);   

  ClusterList newDR({element.GetCluster(badHitType)});
  PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &newDR, "new delta ray", RED);
  PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &newMuonCluster, "NEW MUON ray", VIOLET);
  PandoraMonitoringApi::ViewEvent(this->GetPandora());



  ClusterVector clusterVector;
  clusterVector.push_back(newMuonCluster.front()); clusterVector.push_back(element.GetCluster(badHitType));
  PfoVector pfoVector;
  pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front()); pfoVector.push_back(nullptr);

  pAlgorithm->UpdateForNewCluster(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayRemovalTool::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
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

float DeltaRayRemovalTool::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
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

CartesianVector DeltaRayRemovalTool::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector, const Cluster *const pCluster) const
{
  float shortestDistanceSquared(std::numeric_limits<float>::max());

  CartesianVector closestPoint(0.f,0.f,0.f);
  for (const CartesianVector &testPosition : cartesianPointVector)
  {
    if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > 0.5f)
      continue;

    float separationSquared((testPosition - referencePoint).GetMagnitude());

    if (separationSquared > 10.f)
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

void DeltaRayRemovalTool::ProjectPositions(const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianPointVector &projectedPositions) const
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
    
StatusCode DeltaRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapWindow", m_xOverlapWindow));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

/*

    CaloHitList muonCaloHitList;
    muonClusterList.front()->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);


    CartesianVector deltaRayPoint(0.f,0.f,0.f), muonPoint(0.f,0.f,0.f);
    LArClusterHelper::GetClosestPositions(element.GetCluster(badHitType), muonClusterList.front(), deltaRayPoint, muonPoint);


    // find closest muon projected hits to DR <-- if that ratio is small
    if (projectedHitsFraction < 1.f)
    {
        float closestDistance1(std::numeric_limits<float>::max()), closestDistance2(std::numeric_limits<float>::max());
        CartesianVector closestPoint1(0.f,0.f,0.f), closestPoint2(0.f,0.f,0.f);
        for (const CartesianVector &projectedPosition : muonProjectedPositions)
	{

	  if (LArClusterHelper::GetClosestDistance(projectedPosition, muonClusterList.front()) > 0.5f)
	      continue;

	  float separation((projectedPosition - muonPoint).GetMagnitude());
	  if ((separation > 5.f) && (separation < 20.f))
	    {
	        if(separation < closestDistance1)
		{
		    closestDistance1 = separation;
		    closestPoint1 = projectedPosition;
		}
	    }
	}
        for (const CartesianVector &projectedPosition : muonProjectedPositions)
	{
	    if (LArClusterHelper::GetClosestDistance(projectedPosition, muonClusterList.front()) > 0.5f)
	        continue;

	    float separation((projectedPosition - muonPoint).GetMagnitude());
	    if ((separation > 5.f) && (separation < 20.f))
	    {
		if (separation < closestDistance2)
		{
		    if ((closestPoint1 - projectedPosition).GetMagnitude() > 10.f)
		    {
		        closestDistance2 = separation;
			closestPoint2 = projectedPosition;
		    }
		}
	    }
	}

	if (closestPoint2.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
	{
	  std::cout << "did not found points" << std::endl;
	    return;
	}

        float closestDistance3(std::numeric_limits<float>::max()), closestDistance4(std::numeric_limits<float>::max());
        CartesianVector closestPoint3(0.f,0.f,0.f), closestPoint4(0.f,0.f,0.f);
        for (const CartesianVector &projectedPosition : muonProjectedPositions)
	{
	    if (LArClusterHelper::GetClosestDistance(projectedPosition, muonClusterList.front()) > 0.5f)
	        continue;

	    float separation((closestPoint1 - projectedPosition).GetMagnitude());

	    if (separation < std::numeric_limits<float>::epsilon())
	        continue;

	    if ((projectedPosition - closestPoint2).GetMagnitude() < std::numeric_limits<float>::epsilon())
	        continue;

	    if (separation < closestDistance3)
	    {
	      closestDistance3 = separation;
	      closestPoint3 = projectedPosition;
	    }
	}

        for (const CartesianVector &projectedPosition : muonProjectedPositions)
	{
	    if (LArClusterHelper::GetClosestDistance(projectedPosition, muonClusterList.front()) > 0.5f)
	        continue;

	    float separation((closestPoint2 - projectedPosition).GetMagnitude());

	    if (separation < std::numeric_limits<float>::epsilon())
	        continue;

	    if ((projectedPosition - closestPoint1).GetMagnitude() < std::numeric_limits<float>::epsilon())
	        continue;

	    if (separation < closestDistance4)
	    {
	      closestDistance4 = separation;
	      closestPoint4 = projectedPosition;
	    }
	}

	closestPoint1 = (closestPoint1 + closestPoint3)*(0.5f);
	closestPoint2 = (closestPoint2 + closestPoint4)*(0.5f);

	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint1, "point1", BLACK, 2);
	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPoint2, "point2", BLACK, 2);
	PandoraMonitoringApi::ViewEvent(this->GetPandora());
	std::cout << "found points" << std::endl;

	CartesianVector lineDirection(closestPoint1 - closestPoint2);
	lineDirection = lineDirection.GetUnitVector();
*/
