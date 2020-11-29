/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoViewDeltaRayMatchingAlgorithm::TwoViewDeltaRayMatchingAlgorithm()  :
    m_nMaxMatrixToolRepeats(1000),
    m_minClusterCaloHits(3),
    m_searchRegion1D(3.f),       
    m_xOverlapWindow(1.f),
    m_maxDisplacementSquared(1.0f),
    m_minMatchedFraction(0.5),
    m_minMatchedPoints(2),
    m_pseudoChi2Cut(3.f)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (const Cluster *const pCluster : *pInputClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;
        
        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        selectedClusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::PrepareInputClusters(ClusterList &preparedClusterList)
{
    if (preparedClusterList.empty())
        return;

    const HitType &hitType(LArClusterHelper::GetClusterHitType(preparedClusterList.front()));
    
    this->FillHitToClusterMap(hitType);
    this->FillClusterProximityMap(hitType);
    this->FillClusterToPfoMap(hitType);
    this->FillHitAssociationMap();
}


//------------------------------------------------------------------------------------------------------------------------------------------     
    
void TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const)
{
    TrackTwoViewTopologyOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pCluster1, pCluster2, m_hitAssociationMap, overlapResult));
    
    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------     

StatusCode TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const HitAssociationMap &hitAssociationMap, TrackTwoViewTopologyOverlapResult &overlapResult) const
{
    // Check the overlap of the two clusters
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1); 
    pCluster2->GetClusterSpanX(xMin2, xMax2); 

    const float overlapMinX(std::max(xMin1, xMin2));
    const float overlapMaxX(std::min(xMax1, xMax2));
    const float xOverlap(overlapMaxX - overlapMinX);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    if (!this->AreClustersCompatible(pCluster1, pCluster2))
        return STATUS_CODE_NOT_FOUND;

    // GetProjectedPositions
    CartesianPointVector projectedPositions;
    if (this->GetProjectedPositions(pCluster1, pCluster2, projectedPositions) == STATUS_CODE_NOT_FOUND)
        return STATUS_CODE_NOT_FOUND;

    CaloHitList projectedHits;
    HitOwnershipMap hitOwnershipMap;    
    this->CollectHits(projectedPositions, hitAssociationMap, overlapMinX, overlapMaxX, projectedHits, hitOwnershipMap);
    
    //IGNORE IF OTHER SPANS ARE BETTER <--- They need to also be close to a muon tho (i think it would be if other two are)
    const Cluster *pBestMatchedCluster(nullptr);
    float highestNHits(0);
    for (auto &entry : hitOwnershipMap)
    {
        if (entry.second.size() > highestNHits)
        {
            highestNHits = entry.second.size();
            pBestMatchedCluster = entry.first;
        }
    }
        
    if (pBestMatchedCluster && pBestMatchedCluster->IsAvailable())
    {
        float xMin3(std::numeric_limits<float>::max()), xMax3(-std::numeric_limits<float>::max());
        const OrderedCaloHitList &orderedCaloHitList(pBestMatchedCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                xMin3 = std::min(xMin3, pCaloHit->GetPositionVector().GetX());
                xMax3 = std::max(xMax3, pCaloHit->GetPositionVector().GetX());
            }
        }

        float xOverlap12(xOverlap), xOverlap23(std::min(xMax2, xMax3) - std::max(xMin2, xMin3)), xOverlap13(std::min(xMax1, xMax3) - std::max(xMin1, xMin3));

        if (xOverlap23 > xOverlap12)
            return STATUS_CODE_NOT_FOUND;

        if (xOverlap13 > xOverlap12)
            return STATUS_CODE_NOT_FOUND;
    }

    if (projectedHits.size() == 0)
        return STATUS_CODE_NOT_FOUND;

    /*
    if (this->GetPandora().GetGeometry()->GetLArTPC().GetCenterX() > -370)// && (std::fabs(xOverlap - 2.5313) < 0.001))
    {
        const ClusterList list1({pCluster1}), list2({pCluster2});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &list1, "CLUSTER 1", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &list2, "CLUSTER 2", VIOLET);
    
        for (const CartesianVector &position : projectedPositions)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "projection", RED, 2);

        const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));
        std::cout << "xOverlap: " << xOverlap << std::endl;
        std::cout << "xOverlapFraction: " << xOverlap/xSpan << std::endl;
        std::cout << "nProjectedPoints: " << projectedPositions.size() << std::endl;

        for (const CaloHit *const pCaloHit : projectedHits)
        {
            const CartesianVector &position(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "COLLECTED HIT", DARKGREEN, 2);
        }
        std::cout << "nCollectedHits: " << projectedHits.size() << std::endl;

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    */
    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    ClusterList originalClusterList;
    for (auto &entry : hitOwnershipMap)
        originalClusterList.push_back(entry.first);

    // Create temporary cluster so can calculate overlap result
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_inputClusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));
 
    const Cluster *pTempCluster(nullptr);
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = projectedHits;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pTempCluster));

    // call overlap function
    if (this->CalculateOverlapResult(pCluster1, pCluster2, pTempCluster, projectedHits, overlapResult) == STATUS_CODE_NOT_FOUND)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return STATUS_CODE_NOT_FOUND;
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));

    //overlapResult.SetProjectedCaloHits(projectedHits);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMatchingAlgorithm::AreClustersCompatible(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    ClusterList consideredClusters1, consideredClusters2;
    PfoList nearbyMuonPfos1, nearbyMuonPfos2;
    this->GetNearbyMuonPfos(pCluster1, consideredClusters1, nearbyMuonPfos1);
    this->GetNearbyMuonPfos(pCluster2, consideredClusters2, nearbyMuonPfos2);

    /////////////////////////
    /*
    ClusterList uList({pClusterU}), vList({pClusterV}), wList({pClusterW});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &uList, "pClusterU", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &vList, "pClusterV", DARKGREEN);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &wList, "pClusterW", BLACK);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &nearbyMuonPfosU, "nearbyMuonU", BLUE);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &nearbyMuonPfosV, "nearbyMuonV", ORANGE);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &nearbyMuonPfosW, "nearbyMuonW", VIOLET);
    */
    /////////////////////////
    
    for (const ParticleFlowObject *const pNearbyMuon1 : nearbyMuonPfos1)
    {
        for (const ParticleFlowObject *const pNearbyMuon2 : nearbyMuonPfos2)
        {
            if (pNearbyMuon1 == pNearbyMuon2)
            {
                //std::cout << "FOUND COMMON MUON!" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());                    
                return true;
            }
        }
    }

    //std::cout << "HAVE NOT FOUND COMMON MUON!" << std::endl;
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GetNearbyMuonPfos(const Cluster *const pCluster, ClusterList &consideredClusters, PfoList &nearbyMuonPfos) const 
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const unsigned int hitTypeIndex(m_matchingControl.m_hitTypeToIndexMap.at(hitType));      
    
    const ClusterToPfoMap &clusterToPfoMap((hitTypeIndex == 1) ? m_clusterToPfoMap1 : m_clusterToPfoMap2);
    const ClusterProximityMap &clusterProximityMap((hitTypeIndex == 1) ? m_clusterProximityMap1 : m_clusterProximityMap2);

    consideredClusters.push_back(pCluster);
    
    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pCluster));

    if (clusterProximityIter == clusterProximityMap.end())
        return;
    
    for (const Cluster *const pNearbyCluster : clusterProximityIter->second)
    {
        if (std::find(consideredClusters.begin(), consideredClusters.end(), pNearbyCluster) != consideredClusters.end())
            continue;
        
        const ClusterToPfoMap::const_iterator pfoIter(clusterToPfoMap.find(pNearbyCluster));

        if (pfoIter != clusterToPfoMap.end())
        {
            if (std::find(nearbyMuonPfos.begin(), nearbyMuonPfos.end(), pfoIter->second) == nearbyMuonPfos.end())
                nearbyMuonPfos.push_back(pfoIter->second);
            
            continue;
        }
        
        this->GetNearbyMuonPfos(pNearbyCluster, consideredClusters, nearbyMuonPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::GetProjectedPositions(const Cluster *const pCluster1, const Cluster *const pCluster2,
    CartesianPointVector &projectedPositions) const
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
         return STATUS_CODE_NOT_FOUND;
    
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    if (hitType1 == hitType2)
        return STATUS_CODE_FAILURE;

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
                return statusCodeException.GetStatusCode();

            continue;
        }
    }
    
    const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));
    if ((xOverlap / xSpan < 0.1) || (projectedPositions.size() < 3))
        return STATUS_CODE_NOT_FOUND;
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void TwoViewDeltaRayMatchingAlgorithm::CollectHits(const CartesianPointVector &projectedPositions, const HitAssociationMap &hitAssociationMap,
    const float xMin, const float xMax, CaloHitList &collectedHits, HitOwnershipMap &hitOwnershipMap) const
{
    const ClusterList *pInputCaloHitList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputCaloHitList));

    const PfoList *pDeltaRayList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_deltaRayPfoListName, pDeltaRayList));

    ClusterList deltaRayClusters;
    if(pDeltaRayList)
    {
        for (const ParticleFlowObject *const pPfo : *pDeltaRayList)
            deltaRayClusters.insert(deltaRayClusters.begin(), pPfo->GetClusterList().begin(), pPfo->GetClusterList().end());
    }
    
    for (const CartesianVector &projectedPosition : projectedPositions)
    {
        const CaloHit *pClosestCaloHit(nullptr);
        float closestDistanceSquared(std::numeric_limits<float>::max()), tieBreakerBestEnergy(0.f);
        
        for (const Cluster *const pCluster : *pInputCaloHitList)
        {
            // Protect hits in delta rays
            if (std::find(deltaRayClusters.begin(), deltaRayClusters.end(), pCluster) != deltaRayClusters.end())
                continue;
            
            const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

            for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
            {
                for (const CaloHit *const pCaloHit : *mapEntry.second)
                {
                    const float distanceSquared((pCaloHit->GetPositionVector() - projectedPosition).GetMagnitude());

                    if ((distanceSquared < closestDistanceSquared) ||
                        ((std::fabs(distanceSquared - closestDistanceSquared) < std::numeric_limits<float>::epsilon()) && (pCaloHit->GetHadronicEnergy() > tieBreakerBestEnergy)))
                    {
                        pClosestCaloHit = pCaloHit;
                        closestDistanceSquared = distanceSquared;
                        tieBreakerBestEnergy = pCaloHit->GetHadronicEnergy();
                    }
                }
            }
        }

        if ((closestDistanceSquared < m_maxDisplacementSquared) && (pClosestCaloHit != nullptr) && (std::find(collectedHits.begin(), collectedHits.end(), pClosestCaloHit) == collectedHits.end()))
        {
            this->CollectAssociatedHits(pClosestCaloHit, pClosestCaloHit, hitAssociationMap, xMin, xMax, collectedHits);
            collectedHits.push_back(pClosestCaloHit);
        }
    }

    if (!collectedHits.empty())
    {
        for (const Cluster *const pCluster : *pInputCaloHitList)
        {
            const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

            for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
            {
                for (const CaloHit *const pCaloHit : *mapEntry.second)
                {
                    if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                        hitOwnershipMap[pCluster].push_back(pCaloHit);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------        

void TwoViewDeltaRayMatchingAlgorithm::CollectAssociatedHits(const CaloHit *const pSeedCaloHit, const CaloHit *const pCurrentCaloHit,
    const HitAssociationMap &hitAssociationMap, const float xMin, const float xMax, CaloHitList &associatedHitList) const
{
    HitAssociationMap::const_iterator iter1 = hitAssociationMap.find(pCurrentCaloHit);
    if (iter1 == hitAssociationMap.end())
        return;

    CaloHitVector caloHitVector(iter1->second.begin(), iter1->second.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    for (const CaloHit *const pAssociatedCaloHit : caloHitVector)
    {
        if (pAssociatedCaloHit == pSeedCaloHit)
            continue;

        if ((pAssociatedCaloHit->GetPositionVector().GetX() < xMin) || (pAssociatedCaloHit->GetPositionVector().GetX() > xMax))
            continue;
        
        if (associatedHitList.end() != std::find(associatedHitList.begin(), associatedHitList.end(), pAssociatedCaloHit))
            continue;

        associatedHitList.push_back(pAssociatedCaloHit);

        this->CollectAssociatedHits(pSeedCaloHit, pAssociatedCaloHit, hitAssociationMap, xMin, xMax, associatedHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3,
    const CaloHitList &projectedHits, TrackTwoViewTopologyOverlapResult &overlapResult) const
{
    float xMin1(-std::numeric_limits<float>::max()), xMax1(+std::numeric_limits<float>::max());
    float xMin2(-std::numeric_limits<float>::max()), xMax2(+std::numeric_limits<float>::max());
    float xMin3(-std::numeric_limits<float>::max()), xMax3(+std::numeric_limits<float>::max());

    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);
    pCluster3->GetClusterSpanX(xMin3, xMax3);

    // Need to remove the xPitch from calculations to be consistent with view xSpan calculated in the xOverlapObject
    const float xMinCentre(std::max(xMin1, std::max(xMin2, xMin3)));
    const float xMaxCentre(std::min(xMax1, std::min(xMax2, xMax3)));
    const float xCentreOverlap(xMaxCentre - xMinCentre);

    if (xCentreOverlap < std::numeric_limits<float>::epsilon())
         return STATUS_CODE_NOT_FOUND;
    
    // what is the m_xOverlapWindow? - seems to be like a hit width? (an uncertainty in x)
    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, std::max(xMin2, xMin3)) - xPitch);
    const float xMax(std::min(xMax1, std::min(xMax2, xMax3)) + xPitch);
    const float xOverlap(xMax - xMin);

    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

    if (hitType1 == hitType2 ||  hitType1 == hitType3 || hitType2 == hitType3)
        return STATUS_CODE_FAILURE;

    const unsigned int nPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    // Chi2 calculations
    float pseudoChi2Sum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    
    for (unsigned int n = 0; n < nPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMin1(0.f), zMin2(0.f), zMin3(0.f), zMax1(0.f), zMax2(0.f), zMax3(0.f);
            pCluster1->GetClusterSpanZ(xmin, xmax, zMin1, zMax1);
            pCluster2->GetClusterSpanZ(xmin, xmax, zMin2, zMax2);
            pCluster3->GetClusterSpanZ(xmin, xmax, zMin3, zMax3);

            const float z1(0.5f * (zMin1 + zMax1));
            const float z2(0.5f * (zMin2 + zMax2));
            const float z3(0.5f * (zMin3 + zMax3));

            const float dz1(zMax1 - zMin1);
            const float dz2(zMax2 - zMin2);
            const float dz3(zMax3 - zMin3);
            const float dzPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

            const float zproj1(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, z2, z3));
            const float zproj2(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType3, hitType1, z3, z1));
            const float zproj3(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, z1, z2));

            ++nSamplingPoints;

            const float deltaSquared(((z1 - zproj1) * (z1 - zproj1) + (z2 - zproj2) * (z2 - zproj2) + (z3 - zproj3) * (z3 - zproj3)) / 3.f);
            const float sigmaSquared(dz1 * dz1 + dz2 * dz2 + dz3 * dz3 + dzPitch * dzPitch);
            const float pseudoChi2(deltaSquared / sigmaSquared);

            pseudoChi2Sum += pseudoChi2;
            
            if (pseudoChi2 < m_pseudoChi2Cut)
                ++nMatchedSamplingPoints;
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (statusCodeException.GetStatusCode() != STATUS_CODE_NOT_FOUND)
                return statusCodeException.GetStatusCode();

            continue;
        }
    }

    if (nSamplingPoints == 0)
        return STATUS_CODE_NOT_FOUND;

    const float matchedFraction(static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints));

    /*
    if (this->GetPandora().GetGeometry()->GetLArTPC().GetCenterX() > -370)// && (std::fabs(xOverlap - 2.5313) < 0.001))
    {
        std::cout << "nSamplingPoints: " << nSamplingPoints << std::endl;        
        std::cout << "nMatchedSamplingPoints: " << nMatchedSamplingPoints << std::endl;
        std::cout << "reducedChi2: " << static_cast<float>(pseudoChi2Sum)/static_cast<float>(nSamplingPoints) << std::endl;
    }
    */
    
    if ((matchedFraction < m_minMatchedFraction) || (nMatchedSamplingPoints < m_minMatchedPoints))
        return STATUS_CODE_NOT_FOUND;

    TwoViewXOverlap xOverlapObject(xMin1, xMax1, xMin2, xMax2);
    overlapResult = TrackTwoViewTopologyOverlapResult(nMatchedSamplingPoints, nSamplingPoints, pseudoChi2Sum, xOverlapObject, projectedHits);
    
    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (MatrixToolVector::const_iterator toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end(); )
    {
        if ((*toolIter)->Run(this, this->GetMatchingControl().GetOverlapMatrix()))
        {
            toolIter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxMatrixToolRepeats)
                break;
        }
        else
        {
            ++toolIter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::TidyUp()
{
    m_hitToClusterMap1.clear();
    m_hitToClusterMap2.clear();

    m_kdTree1.clear();
    m_kdTree2.clear();

    m_clusterProximityMap1.clear();
    m_clusterProximityMap2.clear();
    
    m_clusterToPfoMap1.clear();
    m_clusterToPfoMap2.clear();
    
    return BaseAlgorithm::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName)); 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DeltaRayPfoListName", m_deltaRayPfoListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "DeltaRayTools", algorithmToolVector));
    
    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DeltaRayMatrixTool *const pDeltaRayMatrixTool(dynamic_cast<DeltaRayMatrixTool*>(*iter));

        if (!pDeltaRayMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayMatrixTool);
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));    
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SearchRegion1D", m_searchRegion1D));      

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverlapWindow", m_xOverlapWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDisplacementSquared", m_maxDisplacementSquared));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedPoints", m_minMatchedPoints));


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));    

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------     

void TwoViewDeltaRayMatchingAlgorithm::FillHitToClusterMap(const HitType &hitType)
{
    const unsigned int hitTypeIndex(m_matchingControl.m_hitTypeToIndexMap.at(hitType));
    HitToClusterMap &hitToClusterMap((hitTypeIndex == 1) ? m_hitToClusterMap1 : m_hitToClusterMap2);

    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    for (const Cluster *const pCluster : inputClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
                hitToClusterMap[pCaloHit] = pCluster;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::FillClusterProximityMap(const HitType &hitType)
{
    const unsigned int hitTypeIndex(m_matchingControl.m_hitTypeToIndexMap.at(hitType));    
    const HitToClusterMap &hitToClusterMap((hitTypeIndex == 1) ? m_hitToClusterMap1 : m_hitToClusterMap2);

    CaloHitList allCaloHits;
    for (auto &entry : hitToClusterMap)
        allCaloHits.push_back(entry.first);
    
    HitKDTree2D &kdTree((hitTypeIndex == 1) ? m_kdTree1 : m_kdTree2);
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    ClusterProximityMap &clusterProximityMap((hitTypeIndex == 1) ? m_clusterProximityMap1 : m_clusterProximityMap2);
    
    for (const Cluster *const pCluster : inputClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));

                HitKDNode2DList found;
                kdTree.search(searchRegionHits, found);

                for (const auto &hit : found)
                {
                    ClusterList  &nearbyClusterList(clusterProximityMap[pCluster]);
                    const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));

                    if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
                        nearbyClusterList.push_back(pNearbyCluster);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::FillClusterToPfoMap(const HitType &hitType)
{
    const unsigned int hitTypeIndex(m_matchingControl.m_hitTypeToIndexMap.at(hitType));      
    ClusterToPfoMap &clusterToPfoMap((hitTypeIndex == 1) ? m_clusterToPfoMap1 : m_clusterToPfoMap2);
    
    const PfoList *pMuonPfoList(nullptr);
    
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_muonPfoListName, pMuonPfoList));

    if ((!pMuonPfoList) || pMuonPfoList->empty())
        return;

    for (const ParticleFlowObject *const pPfo : *pMuonPfoList)
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pPfo, hitType, pfoClusters);
        for (const Cluster *const pCluster : pfoClusters)
            clusterToPfoMap[pCluster] = pPfo;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::FillHitAssociationMap()
{
    // Third view clusters
    const ClusterList *pInputClusterList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));

    const PfoList *pDeltaRayList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_deltaRayPfoListName, pDeltaRayList));

    ClusterList deltaRayClusters;
    if(pDeltaRayList)
    {
        //can make better by considering hit type... but you can make this better anyway        
        for (const ParticleFlowObject *const pPfo : *pDeltaRayList)
            deltaRayClusters.insert(deltaRayClusters.begin(), pPfo->GetClusterList().begin(), pPfo->GetClusterList().end()); 
    }    
    
    CaloHitList availableCaloHitList;
    for (const Cluster *const pCluster : *pInputClusterList)
    {
        // Protect hits in delta rays
        if (std::find(deltaRayClusters.begin(), deltaRayClusters.end(), pCluster) != deltaRayClusters.end())
            continue;
        
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
                availableCaloHitList.push_back(pCaloHit);
        }
    }
    
    for (const CaloHit *const pCaloHitI : availableCaloHitList)
    {
        for (const CaloHit *const pCaloHitJ : availableCaloHitList)
        {
            if (pCaloHitI == pCaloHitJ)
                continue;

            if ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_maxDisplacementSquared)
            {
                CaloHitList &caloHitListI(m_hitAssociationMap[pCaloHitI]);

                if (caloHitListI.end() == std::find(caloHitListI.begin(), caloHitListI.end(), pCaloHitJ))
                    caloHitListI.push_back(pCaloHitJ);

                CaloHitList &caloHitListJ(m_hitAssociationMap[pCaloHitI]);

                if (caloHitListJ.end() == std::find(caloHitListJ.begin(), caloHitListJ.end(), pCaloHitI))
                    caloHitListJ.push_back(pCaloHitI);

            }
        }
    }
}
    

} // namespace lar_content
