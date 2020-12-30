/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.cc
 *
 *  @brief  Implementation of the delta ray removal tool class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

DeltaRayRemovalTool::DeltaRayRemovalTool() :
    RemovalBaseTool::RemovalBaseTool(),
    m_minSeparation(1.f)
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

        this->RemoveDeltaRayHits(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::RemoveDeltaRayHits(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
    bool &changesMade) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    ClusterSet modifiedClusters, checkedClusters;
    
    for (const TensorType::Element &element : elementList)
    {
        for (const HitType &hitType : hitTypeVector)
	    {
            const Cluster *const pDeltaRayCluster(element.GetCluster(hitType));
            
            if (checkedClusters.count(pDeltaRayCluster))
                continue;
            
            if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) || (modifiedClusters.count(element.GetClusterW())))
                continue;

            if (!this->PassElementChecks(element, hitType))
                continue;

            if (!this->IsContaminated(element, hitType))
                continue;            

            if (!this->IsBestElement(element, hitType, elementList))
                continue;
            
            checkedClusters.insert(pDeltaRayCluster);

	        // Attempt to pull delta ray hits out of muon cluster
            CaloHitList deltaRayHits;
            if (this->CollectDeltaRayHits(element, hitType, deltaRayHits) != STATUS_CODE_SUCCESS)
                continue;
            
            if (deltaRayHits.empty())
                continue;

            modifiedClusters.insert(pDeltaRayCluster);

            this->SplitCluster(pAlgorithm, element, hitType, deltaRayHits);
        }
    }

    changesMade = modifiedClusters.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::PassElementChecks(const TensorType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    if (this->IsMuonEndpoint(element))
        return false;

    const float separation(LArClusterHelper::GetClosestDistance(pDeltaRayCluster, pMuonCluster));

    if (separation > m_minSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::IsMuonEndpoint(const TensorType::Element &element) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    for (const HitType &hitType : hitTypeVector)
    {
        const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
        if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
            return true;

        float xMinDR(+std::numeric_limits<float>::max()), xMaxDR(-std::numeric_limits<float>::max());
        float xMinCR(-std::numeric_limits<float>::max()), xMaxCR(+std::numeric_limits<float>::max());
        
        pDeltaRayCluster->GetClusterSpanX(xMinDR, xMaxDR);
        pMuonCluster->GetClusterSpanX(xMinCR, xMaxCR);

        if ((xMinDR < xMinCR) || (xMaxDR > xMaxCR))
            return true;
	
        float zMinDR(+std::numeric_limits<float>::max()), zMaxDR(-std::numeric_limits<float>::max());
        float zMinCR(-std::numeric_limits<float>::max()), zMaxCR(+std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanZ(xMinDR, xMaxDR, zMinDR, zMaxDR);
        pMuonCluster->GetClusterSpanZ(xMinCR, xMaxCR, zMinCR, zMaxCR);

        if ((zMinDR < zMinCR) || (zMaxDR > zMaxCR))
            return true;
    }

    return false;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::IsContaminated(const TensorType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 100000, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const bool isTransverse(((muonDirection.GetOpeningAngle(xAxis) * 180.f / 3.14) < 20.f) || (((180 - muonDirection.GetOpeningAngle(xAxis)) * 180.f / 3.14) < 20.f));
    
    if (isTransverse)
        return false;

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);

    CaloHitList minusMuonHits, minusDeltaRayHits, plusMuonHits, plusDeltaRayHits;    
    const CartesianVector minusPosition(muonVertex - (muonDirection * 5.f)), plusPosition(muonVertex + (muonDirection * 5.f));

    this->FindExtrapolatedHits(pMuonCluster, muonVertex, minusPosition, minusMuonHits);
    this->FindExtrapolatedHits(pMuonCluster, muonVertex, plusPosition, plusMuonHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, minusPosition, minusDeltaRayHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, plusPosition, plusDeltaRayHits);
    
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

StatusCode DeltaRayRemovalTool::CollectDeltaRayHits(const TensorType::Element &element, const HitType &hitType, CaloHitList &collectedHits) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;
    
    CartesianPointVector muonProjectedPositions;
    this->ProjectMuonPositions(element, hitType, muonProjectedPositions);    

    const float projectedHitsFraction(static_cast<float>(muonProjectedPositions.size()) / pMuonCluster->GetNCaloHits());

    CartesianVector muonDirection(0.f, 0.f, 0.f), positionOnMuon(0.f, 0.f, 0.f);
    if (projectedHitsFraction < 0.8f)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 40, slidingFitPitch);

        CartesianVector deltaRayVertex(0.f,0.f,0.f), muonVertex(0.f,0.f,0.f);
        LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);

        positionOnMuon = this->GetClosestPosition(muonVertex, muonProjectedPositions, pMuonCluster);

	    if (positionOnMuon.GetMagnitude() < std::numeric_limits<float>::epsilon())
            return STATUS_CODE_NOT_FOUND;

        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(positionOnMuon, rL, rT);
        slidingFitResult.GetGlobalFitDirection(rL, muonDirection);
    }

    CaloHitList muonCaloHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);

    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

        for (const CaloHit *const pCaloHit : muonCaloHitList)
        {
            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            
            const float distanceToDeltaRayHits(std::min(LArClusterHelper::GetClosestDistance(hitPosition, pDeltaRayCluster),
                this->GetClosestDistance(pCaloHit, collectedHits)));
            const float distanceToMuonHits((projectedHitsFraction < 0.8f) ? muonDirection.GetCrossProduct(hitPosition - positionOnMuon).GetMagnitude() :
                this->GetClosestDistance(pCaloHit, muonProjectedPositions));

	        if ((std::fabs(distanceToMuonHits - distanceToDeltaRayHits) > std::numeric_limits<float>::epsilon()) && (distanceToDeltaRayHits < distanceToMuonHits)
                    && (distanceToMuonHits > 1.f) && (distanceToDeltaRayHits < 1.f))
	        {
                collectedHits.push_back(pCaloHit);
                hitsAdded = true;
            }
        }
    }

    // Catch if delta ray has travelled along muon
    if ((static_cast<float>(collectedHits.size()) / muonCaloHitList.size()) > 0.05)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const HitType &hitType, CaloHitList &collectedHits) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return;

    pAlgorithm->UpdateUponDeletion(pMuonCluster);    
    pAlgorithm->UpdateUponDeletion(pDeltaRayCluster);
  
    CaloHitList muonCaloHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);
    
    for (const CaloHit *const pCaloHit : muonCaloHitList)
    {
        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, pMuonCluster, pCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*pAlgorithm, pDeltaRayCluster, pCaloHit));
        }
    }

    ClusterVector clusterVector; PfoVector pfoVector;
    clusterVector.push_back(pMuonCluster); clusterVector.push_back(pDeltaRayCluster);
    pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front()); pfoVector.push_back(nullptr);

    pAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode DeltaRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

