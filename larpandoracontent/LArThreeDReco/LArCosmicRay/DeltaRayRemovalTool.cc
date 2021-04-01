/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.cc
 *
 *  @brief  Implementation of the delta ray removal tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.h"

using namespace pandora;

namespace lar_content
{

DeltaRayRemovalTool::DeltaRayRemovalTool() :
    m_minSeparation(1.f),
    m_slidingFitWindow(10000),
    m_minDeviationFromTransverse(0.35f),
    m_contaminationWindow(5.f),
    m_significantHitThreshold(3),
    m_minDistanceFromMuon(1.f),
    m_maxDistanceToCollected(1.f)
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

        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList);
        
        for (const TensorType::Element &element : elementList)
            usedKeyClusters.insert(element.GetClusterU());

        this->ExamineConnectedElements(pAlgorithm, elementList, changesMade);
    }
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
    bool &changesMade) const
{
    ClusterSet modifiedClusters, checkedClusters;
    
    for (const TensorType::Element &element : elementList)
    {
        for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
	    {
            const Cluster *pDeltaRayCluster(element.GetCluster(hitType));
            const ParticleFlowObject *const pMuonPfo(element.GetOverlapResult().GetCommonMuonPfoList().front());

            if (checkedClusters.count(pDeltaRayCluster))
                continue;
            
            if ((modifiedClusters.count(element.GetClusterU())) || (modifiedClusters.count(element.GetClusterV())) || (modifiedClusters.count(element.GetClusterW())))
                continue;

            if (!this->PassElementChecks(pAlgorithm, element, hitType))
                continue;

            if (!this->IsContaminated(pAlgorithm, element, hitType))
                continue;            

            if (!this->IsBestElement(pAlgorithm, element, hitType, elementList))
                continue;
            
            checkedClusters.insert(pDeltaRayCluster);

            CaloHitList deltaRayHits;
            if (pAlgorithm->CollectHitsFromMuon(nullptr, nullptr, pDeltaRayCluster, pMuonPfo, m_minDistanceFromMuon, m_maxDistanceToCollected, deltaRayHits) != STATUS_CODE_SUCCESS)
                continue;
            
            modifiedClusters.insert(pDeltaRayCluster);

            this->SplitMuonCluster(pAlgorithm, element, hitType, deltaRayHits);
         }
    }

    changesMade = modifiedClusters.size();
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::PassElementChecks(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &hitType) const
{
    // ATTN: Avoid endpoints, topological michel reconstruction is very ambiguous
    if (this->IsMuonEndpoint(pAlgorithm, element, false))
        return false;

    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;    

    const float separation(LArClusterHelper::GetClosestDistance(pDeltaRayCluster, pMuonCluster));

    if (separation > m_minSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRemovalTool::IsContaminated(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &hitType) const
{
    const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        return false;
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, m_slidingFitWindow, slidingFitPitch);

    CartesianVector muonDirection(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalDirection(slidingFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), muonDirection);

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const bool isTransverse((muonDirection.GetOpeningAngle(xAxis) < m_minDeviationFromTransverse) ||
        ((M_PI - muonDirection.GetOpeningAngle(xAxis)) < m_minDeviationFromTransverse));
    
    if (isTransverse)
        return false;

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArClusterHelper::GetClosestPositions(pDeltaRayCluster, pMuonCluster, deltaRayVertex, muonVertex);

    CaloHitList minusMuonHits, minusDeltaRayHits, plusMuonHits, plusDeltaRayHits;    
    const CartesianVector minusPosition(muonVertex - (muonDirection * m_contaminationWindow)); 
    const CartesianVector plusPosition(muonVertex + (muonDirection * m_contaminationWindow));

    this->FindExtrapolatedHits(pMuonCluster, muonVertex, minusPosition, minusMuonHits);
    this->FindExtrapolatedHits(pMuonCluster, muonVertex, plusPosition, plusMuonHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, minusPosition, minusDeltaRayHits);
    this->FindExtrapolatedHits(pDeltaRayCluster, muonVertex, plusPosition, plusDeltaRayHits);

    if ((minusMuonHits.size() < m_significantHitThreshold) && (plusMuonHits.size() < m_significantHitThreshold))
        return true;

    if ((minusMuonHits.size() < m_significantHitThreshold) && (minusDeltaRayHits.size() >= m_significantHitThreshold) &&
        (plusDeltaRayHits.size() < m_significantHitThreshold) && (plusMuonHits.size() >= m_significantHitThreshold))
    {
        return true;
    }

    if ((plusMuonHits.size() < m_significantHitThreshold) && (plusDeltaRayHits.size() >= m_significantHitThreshold) &&
        (minusDeltaRayHits.size() < m_significantHitThreshold)  && (minusMuonHits.size() >= m_significantHitThreshold))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRemovalTool::SplitMuonCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const HitType &hitType, const CaloHitList &deltaRayHits) const
{
    const Cluster *pDeltaRayCluster(element.GetCluster(hitType)), *pMuonCluster(nullptr);
    
    if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    pAlgorithm->UpdateUponDeletion(pMuonCluster); pAlgorithm->UpdateUponDeletion(pDeltaRayCluster);
    
    pAlgorithm->SplitMuonCluster(pAlgorithm->GetClusterListName(hitType), pMuonCluster, deltaRayHits, pDeltaRayCluster);

    ClusterVector clusterVector; PfoVector pfoVector;
    clusterVector.push_back(pMuonCluster); pfoVector.push_back(element.GetOverlapResult().GetCommonMuonPfoList().front());
    clusterVector.push_back(pDeltaRayCluster); pfoVector.push_back(nullptr);
    pAlgorithm->UpdateForNewClusters(clusterVector, pfoVector); 
}    

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode DeltaRayRemovalTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDeviationFromTransverse", m_minDeviationFromTransverse));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ContaminationWindow", m_contaminationWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SignificantHitThreshold", m_significantHitThreshold));      

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDistanceFromMuon", m_minDistanceFromMuon));      

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceToCollected", m_maxDistanceToCollected));      

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

