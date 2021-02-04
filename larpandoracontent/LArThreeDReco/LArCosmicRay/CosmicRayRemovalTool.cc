/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CosmicRayRemovalTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayRemovalTool.h"

using namespace pandora;

namespace lar_content
{

CosmicRayRemovalTool::CosmicRayRemovalTool() :
    RemovalBaseTool::RemovalBaseTool(),    
    m_minSeparation(2.f),
    m_minContaminationLength(3.f)
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

        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList);
        
        for (const TensorType::Element &element : elementList)
            usedKeyClusters.insert(element.GetCluster(TPC_VIEW_U));

        this->ExamineConnectedElements(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &changesMade) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    ClusterSet modifiedClusters, checkedClusters;
    
    for (const TensorType::Element &element : elementList)
    {        
        for (const HitType &hitType : hitTypeVector)
	    {
            if (checkedClusters.count(element.GetCluster(hitType)))
                continue;
            
            if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) || (modifiedClusters.count(element.GetClusterW())))
                continue;
            
            if (!this->PassElementChecks(element, hitType))
                continue;
            
            if (!this->IsContaminated(element, hitType))
                continue;

            if (!this->IsBestElement(element, hitType, elementList))
                continue;

            // ATTN: For speed optimisation
            checkedClusters.insert(element.GetCluster(hitType));

	        // Attempt to pull delta ray hits out of cluster
            CaloHitList deltaRayHits;
            this->CreateSeed(pAlgorithm, element, hitType, deltaRayHits);

            if (deltaRayHits.empty())
                continue;

            // If seed cannot be grown, abort
            CaloHitList remnantHits;
            if (this->GrowSeed(pAlgorithm, element, hitType, deltaRayHits, remnantHits) != STATUS_CODE_SUCCESS)
                continue;

            // If found all delta ray hits, abort
            if (deltaRayHits.size() == element.GetCluster(hitType)->GetNCaloHits())
                continue;

            modifiedClusters.insert(element.GetCluster(hitType));

            this->SplitCluster(pAlgorithm, element, hitType, deltaRayHits, remnantHits);
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

    if (this->IsMichel(element, hitType))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool CosmicRayRemovalTool::IsMichel(const TensorType::Element &element, const HitType &hitType) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &otherHitType : hitTypeVector)
    {
        if (otherHitType == hitType)
            continue;

        const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
        if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
            return true;

        float xMinDR(-std::numeric_limits<float>::max()), xMaxDR(+std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanX(xMinDR, xMaxDR);
        
        float xMinCR(-std::numeric_limits<float>::max()), xMaxCR(+std::numeric_limits<float>::max());
        pMuonCluster->GetClusterSpanX(xMinCR, xMaxCR);

        if ((xMinDR < xMinCR) || (xMaxDR > xMaxCR))
            return true;
	
        float zMinDR(-std::numeric_limits<float>::max()), zMaxDR(+std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanZ(xMinDR, xMaxDR, zMinDR, zMaxDR);
        
        float zMinCR(-std::numeric_limits<float>::max()), zMaxCR(+std::numeric_limits<float>::max());
        pMuonCluster->GetClusterSpanZ(xMinCR, xMaxCR, zMinCR, zMaxCR);

        if ((zMinDR < zMinCR) || (zMaxDR > zMaxCR))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayRemovalTool::IsContaminated(const TensorType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);    

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 10000, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    // Find furthest point in DR cluster that lies on the projected muon track
    float furthestSeparation(0.f);
    CartesianVector extendedPoint(0.f, 0.f, 0.f);
        
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

//LEAVING THIS HERE BECAUSE I THINK THIS MAY CAUSE ISSUES IN THE SLIDING FITS OF PFOS
void CosmicRayRemovalTool::CreateSeed(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const HitType &hitType, CaloHitList &collectedHits) const
{
    CartesianPointVector muonProjectedPositions;
    if (pAlgorithm->ProjectMuonPositions(hitType, element.GetOverlapResult().GetCommonMuonPfoList().front(), muonProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CartesianPointVector deltaRayProjectedPositions;   
    if (this->ProjectDeltaRayPositions(pAlgorithm, element, hitType, deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
        return;

    CaloHitList deltaRayHitList;
    element.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);
    
    for (const CaloHit *const pCaloHit : deltaRayHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());

        for (const CartesianVector &projectedPosition : deltaRayProjectedPositions)
        {
            const float distanceToProjectionSquared((position - projectedPosition).GetMagnitudeSquared());

            if (distanceToProjectionSquared < 1.f * 1.f)
	        {
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

StatusCode CosmicRayRemovalTool::GrowSeed(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const HitType &hitType, CaloHitList &deltaRayHits, CaloHitList &remnantHits) const
{
    const Cluster *const pDeltaRayCluster(element.GetCluster(hitType)), *pMuonCluster(nullptr);
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    CartesianPointVector muonProjectedPositions;
    if (pAlgorithm->ProjectMuonPositions(hitType, element.GetOverlapResult().GetCommonMuonPfoList().front(), muonProjectedPositions) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    // To avoid fluctuatuions, parameterise the muon track
    CartesianVector muonDirection(0.f, 0.f, 0.f), positionOnMuon(0.f, 0.f, 0.f);
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 40, slidingFitPitch);

    CartesianVector deltaRayVertex(0.f,0.f,0.f), muonVertex(0.f,0.f,0.f);
    LArClusterHelper::GetClosestPositions(element.GetCluster(hitType), pMuonCluster, deltaRayVertex, muonVertex);

    positionOnMuon = LArMuonLeadingHelper::GetClosestPosition(muonVertex, muonProjectedPositions, pMuonCluster);

    if (positionOnMuon.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(positionOnMuon, rL, rT);
    slidingFitResult.GetGlobalFitDirection(rL, muonDirection);

    // Identify delta ray hits
    this->CollectHits(positionOnMuon, muonDirection, pDeltaRayCluster, 0.5f, true, deltaRayHits, deltaRayHits);

    // Identify remnant hits
    this->CollectHits(positionOnMuon, muonDirection, pDeltaRayCluster, 1.f, false, deltaRayHits, remnantHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayRemovalTool::CollectHits(const CartesianVector &positionOnMuon, const CartesianVector &muonDirection, const Cluster *const pDeltaRayCluster,
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

void CosmicRayRemovalTool::SplitCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &hitType,
    CaloHitList &collectedHits, CaloHitList &deltaRayRemnantHits) const
{
    const Cluster *pMuonCluster(nullptr);
    
    if (this->GetMuonCluster(element, hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return;

    pAlgorithm->UpdateUponDeletion(pMuonCluster);    
    pAlgorithm->UpdateUponDeletion(element.GetCluster(hitType));

    std::string originalListName, fragmentListName;
    ClusterList originalClusterList(1, element.GetCluster(hitType));
    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetClusterListName(hitType)));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, originalClusterList, originalListName, fragmentListName));

    CaloHitList deltaRayHitList;
    element.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

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

void CosmicRayRemovalTool::FragmentRemnant(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const HitType &hitType, const Cluster *const pMuonCluster,
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

StatusCode CosmicRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinContaminationLength", m_minContaminationLength));    

    return RemovalBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content

