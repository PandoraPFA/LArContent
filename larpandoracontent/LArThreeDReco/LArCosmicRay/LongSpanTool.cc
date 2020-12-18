/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/LongSpanTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/LongSpanTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

LongSpanTool::LongSpanTool() :
    m_xOverlapWindow(0.1f),
    m_minXOverlapFraction(0.7f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongSpanTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->InvestigateLongSpans(pAlgorithm, elementList, changesMade);
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongSpanTool::InvestigateLongSpans(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
    for (TensorType::Element &element : elementList)
    {
        const Cluster *pLongCluster(nullptr);

        // Check if long span creates ambiguities
        if(!this->GetLongCluster(element, pLongCluster))
            continue;

	const HitType &badHitType(LArClusterHelper::GetClusterHitType(pLongCluster));

        ////////////////////////// 
	/*        
        std::cout << "uSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_U) << std::endl;
        std::cout << "vSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_V) << std::endl;
        std::cout << "wSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_W) << std::endl;
        std::cout << "overlap: " << element.GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;
                
        ClusterList uCluster({element.GetCluster(TPC_VIEW_U)}), vCluster({element.GetCluster(TPC_VIEW_V)}), wCluster({element.GetCluster(TPC_VIEW_W)});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &uCluster, "uCluster_1", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &vCluster, "vCluster_1", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &wCluster, "wCluster_1", VIOLET);

        ClusterList longCluster({pLongCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &longCluster, "longCluster", BLACK);
	*/

        ////////////////////////// 
	// Also checks only one common muon
        if (!this->IsConnected(element))
        {
            std::cout << "reject because not enough connection points" << std::endl;
	    // PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }

        //PandoraMonitoringApi::ViewEvent(this->GetPandora()); 

	// Check if ambiguity is caused by muon issues
	if (this->IsMuonEndpoint(element, badHitType))
	{
            //////////////////////// 

        std::cout << "uSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_U) << std::endl;
        std::cout << "vSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_V) << std::endl;
        std::cout << "wSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_W) << std::endl;
        std::cout << "overlap: " << element.GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;
                
        ClusterList uCluster({element.GetCluster(TPC_VIEW_U)}), vCluster({element.GetCluster(TPC_VIEW_V)}), wCluster({element.GetCluster(TPC_VIEW_W)});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &uCluster, "uCluster_1", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &vCluster, "vCluster_1", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &wCluster, "wCluster_1", VIOLET);

        ClusterList longCluster({pLongCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &longCluster, "longCluster", BLACK);
	PandoraMonitoringApi::ViewEvent(this->GetPandora());
	    ////////////////////////// 

	    // Attempt to pull delta ray hits out of cluster
            CaloHitList collectedHits;
	    this->CreateSeed(element, badHitType, collectedHits);

            ////////////////////////// 


            ////////////////////////// 
	    this->GrowSeed(element, badHitType, collectedHits);

            for (const CaloHit *const pCaloHit : collectedHits)
            {
	      const CartesianVector &position(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "collected", VIOLET, 2);
            }

	    if ((!collectedHits.empty()) && (collectedHits.size() != element.GetCluster(badHitType)->GetNCaloHits()))
	    {
	      this->SplitCluster(pAlgorithm, element, badHitType, collectedHits);
	      changesMade = true;
	    }
	}
	else
	{
   	    // asumme all is well and investigate other views for missing stray clusters
	    float spanMinX(element.GetOverlapResult().GetViewMinX(badHitType));
	    float spanMaxX(element.GetOverlapResult().GetViewMaxX(badHitType));

	    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
	    for (const HitType &hitType : hitTypeVector)
	    {
	        if (hitType == badHitType)
		  continue;

		ClusterList collectedClusters;
		pAlgorithm->CollectStrayHits(element.GetCluster(hitType), spanMinX, spanMaxX, collectedClusters);

                // Add stray clusters in
		if (!collectedClusters.empty())
                {
		  std::cout << "HERE" << std::endl;
		    pAlgorithm->AddInStrayClusters(element.GetCluster(hitType), collectedClusters);
		    changesMade = true;
		}
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongSpanTool::GetLongCluster(const TensorType::Element &element, const Cluster *&pLongCluster) const
{
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    float shortestSpan(std::numeric_limits<float>::max()), longestSpan(-std::numeric_limits<float>::max());
    for (const HitType &hitType : hitTypeVector)
    {
        const float viewXSpan(element.GetOverlapResult().GetViewXSpan(hitType));
        
        if (viewXSpan > longestSpan)
        {
            longestSpan = viewXSpan;
            pLongCluster = element.GetCluster(hitType);
        }

        if (viewXSpan < shortestSpan)
            shortestSpan = viewXSpan;
    }

    // Check if failed because of span inconsistencies
    if ((longestSpan < std::numeric_limits<float>::epsilon()) || (element.GetOverlapResult().GetXOverlap().GetXOverlapSpan() / longestSpan > m_minXOverlapFraction))
    {
        std::cout << "DID NOT FAIL BECUASE OF SPAN" << std::endl;
        return false;
    }

    float middleSpan(std::numeric_limits<float>::max());
    for (const HitType &hitType : hitTypeVector)
    {
        const float viewXSpan(element.GetOverlapResult().GetViewXSpan(hitType));

        if (std::fabs(viewXSpan - shortestSpan) < std::numeric_limits<float>::epsilon())
            continue;

        if (std::fabs(viewXSpan - longestSpan) < std::numeric_limits<float>::epsilon())
            continue;
        
        middleSpan = viewXSpan;
    }

    // Was the failure a result of a long span?
    return ((longestSpan - middleSpan) > (middleSpan - shortestSpan));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongSpanTool::IsConnected(const TensorType::Element &element) const
{
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
        return false;

    unsigned int connectedClusterCount(0);
    for (const HitType &hitType : hitTypeVector)
    {
        ClusterList muonClusterList;
        LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);

        if (muonClusterList.size() != 1)
        {
            std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
            continue;
        }

        const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList));

        if (separation < 2.f)
            ++connectedClusterCount;
    }

    return (connectedClusterCount > 1);

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongSpanTool::IsMuonEndpoint(const TensorType::Element &element, const HitType &badHitType) const
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), badHitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
    {
        std::cout << "ISOBEL MUON CLUSTER SIZE DOES NOT EQUAL ONE" << std::endl;
        return false;
    }

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(muonClusterList.front(), 20, slidingFitPitch);

    const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(badHitType), muonClusterList));
    
    // If delta ray has a connection point on the muon
    if (separation < 1.f)
    {
        CartesianVector deltaRayPoint(0.f,0.f,0.f), muonPoint(0.f,0.f,0.f);
        LArClusterHelper::GetClosestPositions(element.GetCluster(badHitType), muonClusterList.front(), deltaRayPoint, muonPoint);

        return (this->ShouldSplitDeltaRay(muonClusterList.front(), element.GetCluster(badHitType), muonPoint, slidingFitResult));
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongSpanTool::SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
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

  /*
  std::string originalListName, fragmentListName;
  ClusterList originalClusterList(1, element.GetCluster(badHitType));

  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(badHitType)));
  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, originalClusterList, originalListName, fragmentListName));
 
  const Cluster *pDeltaRay(nullptr);
  PandoraContentApi::Cluster::Parameters parameters;
  parameters.m_caloHitList = collectedHits;
  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, parameters, pDeltaRay));

  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*pAlgorithm, fragmentListName, originalListName));
  */

  ClusterList newMuonCluster;
  LArPfoHelper::GetClusters(element.GetOverlapResult().GetCommonMuonPfoList().front(), badHitType, newMuonCluster);   

  ClusterVector clusterVector;
  clusterVector.push_back(newMuonCluster.front()); clusterVector.push_back(element.GetCluster(badHitType));//clusterVector.push_back(pDeltaRay);
  PfoVector pfoVector;
  pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front()); pfoVector.push_back(nullptr);

  pAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongSpanTool::GrowSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

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
	      return;
	    }

	    break;
        }
	break;
    }
            
    CartesianPointVector muonProjectedPositions;
    this->ProjectPositions(muonClusterList1.front(), muonClusterList2.front(), muonProjectedPositions);

    for (const CartesianVector &position : muonProjectedPositions)
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", BLACK, 2);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    CaloHitList longClusterCaloHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(longClusterCaloHitList);

    bool hitsAdded(true);
    while (hitsAdded)
    {
      hitsAdded = false;

        for (const CaloHit *const pCaloHit : longClusterCaloHitList)
        {
	    if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
	        continue;

	    float distanceToDeltaRayHitsSquared(this->GetClosestDistance(pCaloHit, collectedHits));
	    float distanceToMuonHitsSquared(this->GetClosestDistance(pCaloHit, muonProjectedPositions));

	    if ((std::fabs(distanceToMuonHitsSquared - distanceToDeltaRayHitsSquared) > std::numeric_limits<float>::epsilon()) && (distanceToDeltaRayHitsSquared < distanceToMuonHitsSquared)
		&& distanceToMuonHitsSquared > 1.0f)
	    {
	      collectedHits.push_back(pCaloHit);
	      hitsAdded = true;
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LongSpanTool::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
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

float LongSpanTool::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
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

void LongSpanTool::CreateSeed(const TensorType::Element &element, const HitType &badHitType, CaloHitList &collectedHits) const
{
    CartesianPointVector projectedPositions;
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == badHitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == badHitType) || (hitType1 == hitType2))
                continue;

	    this->ProjectPositions(element.GetCluster(hitType1), element.GetCluster(hitType2), projectedPositions);

	    break;
        }
	break;
    }

    CaloHitList longClusterCaloHitList;
    element.GetCluster(badHitType)->GetOrderedCaloHitList().FillCaloHitList(longClusterCaloHitList);
    
    for (const CaloHit *const pCaloHit : longClusterCaloHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());

        for (const CartesianVector &projectedPosition : projectedPositions)
        {
            const float distanceSquared((position - projectedPosition).GetMagnitude());

            if (distanceSquared < 4.f)
	    {
                collectedHits.push_back(pCaloHit);
		break;
	    }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongSpanTool::ProjectPositions(const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianPointVector &projectedPositions) const
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

bool LongSpanTool::ShouldSplitDeltaRay(const Cluster *const pMuonCluster, const Cluster *const pDeltaRayCluster, const CartesianVector &muonPosition,
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

    if (minusMuonHits.empty() && plusMuonHits.empty())
        return false;

    if (minusDeltaRayHits.empty() && plusDeltaRayHits.empty())
        return false;

    std::cout << "minusMuonHits.size(): " << minusMuonHits.size() << std::endl;
    std::cout << "plusMuonHits.size(): " << plusMuonHits.size() << std::endl;
    std::cout << "minusDeltaRayHits.size(): " << minusDeltaRayHits.size() << std::endl;
    std::cout << "plusDeltaRayHits.size(): " << plusDeltaRayHits.size() << std::endl;

    // change this to be like, is the delta ray hits on the rhs close to that of the muons on the lhs?
    if ((minusMuonHits.size() < 2) && (plusDeltaRayHits.size() < 2) && (minusDeltaRayHits.size() > 2) && (plusMuonHits.size() > 2))
        return true;

    if ((plusMuonHits.size() < 2) && (minusDeltaRayHits.size() < 2) && (plusDeltaRayHits.size() > 2) && (minusMuonHits.size() > 2))
        return true;

    return false;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongSpanTool::FindExtrapolatedHits(const Cluster *const pCluster, const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary,
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

bool LongSpanTool::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
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

bool LongSpanTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, const CartesianVector &lineEnd, const float distanceToLine) const
{
    CartesianVector lineDirection(lineStart - lineEnd);
    lineDirection = lineDirection.GetUnitVector();
    
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
       return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongSpanTool::GetGoodXOverlapExtrema(const TensorType::Element &element, const HitType &badHitType, float &minX, float &maxX) const
{
    if (badHitType == TPC_VIEW_U)
    {
        minX = std::max(element.GetOverlapResult().GetXOverlap().GetVMinX(), element.GetOverlapResult().GetXOverlap().GetWMinX());
        maxX = std::min(element.GetOverlapResult().GetXOverlap().GetVMaxX(), element.GetOverlapResult().GetXOverlap().GetWMaxX());
    }

    if (badHitType == TPC_VIEW_V)
    {
        minX = std::max(element.GetOverlapResult().GetXOverlap().GetUMinX(), element.GetOverlapResult().GetXOverlap().GetWMinX());
        maxX = std::min(element.GetOverlapResult().GetXOverlap().GetUMaxX(), element.GetOverlapResult().GetXOverlap().GetWMaxX());
    }

    if (badHitType == TPC_VIEW_W)
    {
        minX = std::max(element.GetOverlapResult().GetXOverlap().GetUMinX(), element.GetOverlapResult().GetXOverlap().GetVMinX());
        maxX = std::min(element.GetOverlapResult().GetXOverlap().GetUMaxX(), element.GetOverlapResult().GetXOverlap().GetVMaxX());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode LongSpanTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content


/*

  HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
        return; //false

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == badHitType)
            continue;

        for (const HitType &hitType2 : hitTypeVector)
        {
            if (hitType2 == badHitType)
                continue;

            if (hitType1 == hitType2)
                continue;
        
            ClusterList muonClusterList1, muonClusterList2;
            LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType1, muonClusterList1);
            LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType2, muonClusterList2);
            
            if ((muonClusterList1.size() != 1) || (muonClusterList2.size() != 1))
            {
                std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
                continue;
            }
            
            const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            
            const TwoDSlidingFitResult slidingFitResult1(muonClusterList1.front(), 20, slidingFitPitch), slidingFitResult2(muonClusterList2.front(), 20, slidingFitPitch);
            LArPointingCluster pointingCluster1(muonClusterList1.front(), 20, slidingFitPitch), pointingCluster2(muonClusterList2.front(), 20, slidingFitPitch);

            const CartesianVector &innerVertex1(pointingCluster1.GetInnerVertex().GetPosition());
            CartesianVector projectedInnerVertex1(0.f, 0.f, 0.f);
            
            StatusCode innerStatusCode(slidingFitResult2.GetGlobalFitPositionAtX(innerVertex1.GetX(), projectedInnerVertex1));
            if (innerStatusCode == STATUS_CODE_NOT_FOUND)
                slidingFitResult2.GetExtrapolatedPositionAtX(innerVertex1.GetX(), projectedInnerVertex1);

            const CartesianVector &outerVertex1(pointingCluster1.GetOuterVertex().GetPosition());
            CartesianVector projectedOuterVertex1(0.f, 0.f, 0.f);
            
            StatusCode outerStatusCode(slidingFitResult2.GetGlobalFitPositionAtX(outerVertex1.GetX(), projectedOuterVertex1));
            if (outerStatusCode == STATUS_CODE_NOT_FOUND)
                slidingFitResult2.GetExtrapolatedPositionAtX(outerVertex1.GetX(), projectedOuterVertex1);

            // can still be zerro so just ignore that case if it happens?
            
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerVertex1, "innerVertex1", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectedInnerVertex1, "projectedInnerVertex1", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerVertex1, "outerVertex1", BLUE, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectedOuterVertex1, "projectedOuterVertex1", BLUE, 2);            

            std::cout << "innerVertex1: " << innerVertex1 << std::endl;            
            std::cout << "projectedInnerVertex1: " << projectedInnerVertex1 << std::endl;
            std::cout << "outerVertex1: " << outerVertex1 << std::endl;            
            std::cout << "projectedOuterVertex1: " << projectedOuterVertex1 << std::endl;            
            
            PandoraMonitoringApi::Pause(this->GetPandora());

            CartesianVector innerProjection(0.f,0.f,0.f), outerProjection(0.f,0.f,0.f);
            float chi2;
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, innerVertex1, projectedInnerVertex1, innerProjection, chi2);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, outerVertex1, projectedOuterVertex1, outerProjection, chi2);
            
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerProjection, "innerProjection", RED, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerProjection, "outerProjection", VIOLET, 2);
            
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
*/



/*
void TrackExtensionRefinementAlgorithm::GetExtrapolatedCaloHits(const ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList,
    const ClusterList &createdMainTrackClusters, ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    // Construct initial fit
    CartesianVector extrapolatedStartPosition(clusterAssociation.IsEndUpstream() ? downstreamPoint : upstreamPoint);
    CartesianVector extrapolatedDirection(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
    const CartesianVector clusterSubsetBoundary(extrapolatedStartPosition + (extrapolatedDirection * (-1.f) * m_growingFitInitialLength));

    // Collect extrapolated hits by performing a running fit
    unsigned int count(0);
    CartesianVector extrapolatedEndPosition(0.f, 0.f, 0.f);
    unsigned int hitsCollected(std::numeric_limits<int>::max());

    while (hitsCollected)
    {
        hitsCollected = 0;

        extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

        for (const CaloHit *const pCaloHit : hitsInRegion.at(pCluster))
        {
                    // ATTN: To avoid counting same hit twice
                    if (!this->IsInLineSegment(extrapolatedStartPosition, extrapolatedEndPosition, hitPosition))
                        continue;

                    if (!this->IsCloseToLine(hitPosition, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
                        continue;

                    ++hitsCollected;
                    
                    runningFitPositionVector.push_back(hitPosition);
                }
            }
        }
        catch (const StatusCodeException &)
        {
            return;
        }

        ++count;
    }
}
*/

/*

    LArPointingCluster pointingCluster(muonClusterList.front(), 20, slidingFitPitch);

    const CartesianVector &innerVertex(pointingCluster.GetInnerVertex().GetPosition()), &outerVertex(pointingCluster.GetOuterVertex().GetPosition());

    const float innerSeparation(LArClusterHelper::GetClosestDistance(innerVertex, element.GetCluster(badHitType)));
    const float outerSeparation(LArClusterHelper::GetClosestDistance(outerVertex, element.GetCluster(badHitType)));

    if ((innerSeparation < 2.f) || (outerSeparation < 2.f))
    {
        const CartesianVector innerPoint(LArClusterHelper::GetClosestPosition(innerVertex, element.GetCluster(badHitType)));
        const CartesianVector outerPoint(LArClusterHelper::GetClosestPosition(outerVertex, element.GetCluster(badHitType)));

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerPoint, "inner", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerPoint, "OUTER", BLACK, 2);
        
        std::cout << "DELTA RAY IS CLOSE TO MUON ENDPOINT" << std::endl;

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
*/
