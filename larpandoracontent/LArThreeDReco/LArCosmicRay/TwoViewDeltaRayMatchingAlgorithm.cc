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
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pCluster1, pCluster2, overlapResult));
    
    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------     
    
void TwoViewDeltaRayMatchingAlgorithm::RemoveThirdViewCluster(const Cluster *const pCluster)
{
    //std::cout << "REMOVING CLUSTER..." << std::endl;

    //std::cout << "REMOVED CLUSTER: " << pCluster << std::endl;
    
    auto &theMatrix(this->GetMatchingControl().GetOverlapMatrix());

    for (auto iter1 = theMatrix.begin(); iter1 != theMatrix.end(); ++iter1)
    {
        auto iter2(iter1->second);
        
        for (auto &entry : iter2)
        {
            TrackTwoViewTopologyOverlapResult &overlapResult(entry.second);

            ClusterList matchedClusters(overlapResult.GetMatchedClusterList());

            auto matchedClustersIter(std::find(matchedClusters.begin(), matchedClusters.end(), pCluster));

            if (matchedClustersIter == matchedClusters.end())
                continue;
            
            matchedClusters.erase(matchedClustersIter);

            //std::cout << "pCluster1: " << (iter1->first) << std::endl;
            //std::cout << "pCluster2: " << (entry.first) << std::endl;            

            //std::cout << "OVERLAP LIST SIZE BEFORE: " << overlapResult.GetMatchedClusterList().size() << std::endl;
            //overlapResult.SetMatchedClusterList(matchedClusters);
            //std::cout << "OVERLAP LIST SIZE AFTER: " << overlapResult.GetMatchedClusterList().size() << std::endl;            

            const Cluster *pBestMatchedCluster(nullptr); float reducedChiSquared(1.5f);
            this->GetBestMatchedCluster(iter1->first, entry.first, overlapResult.GetCommonMuonPfoList(), matchedClusters, pBestMatchedCluster, reducedChiSquared);

            overlapResult = TrackTwoViewTopologyOverlapResult(overlapResult.GetXOverlap(), overlapResult.GetCommonMuonPfoList(), pBestMatchedCluster, matchedClusters, reducedChiSquared);
            //overlapResult.SetBestMatchedCluster(pBestMatchedCluster);
        }
    }

    //std::cout << "ET FIN" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------     

StatusCode TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, TrackTwoViewTopologyOverlapResult &overlapResult) const
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

    PfoList commonMuonPfoList;
    this->AreClustersCompatible(pCluster1, pCluster2, commonMuonPfoList);
    
    if (commonMuonPfoList.empty())
        return STATUS_CODE_NOT_FOUND;

    // GetProjectedPositions
    CartesianPointVector projectedPositions;
    if (this->GetProjectedPositions(pCluster1, pCluster2, projectedPositions) == STATUS_CODE_NOT_FOUND)
        return STATUS_CODE_NOT_FOUND;

    // ATTN: Find all matched clusters, including unavailable
    ClusterList matchedClusterList;
    this->CollectHits(pCluster1, pCluster2, projectedPositions, matchedClusterList);

    if (matchedClusterList.empty())
        return STATUS_CODE_NOT_FOUND;

    float reducedChiSquared(1.5f);    
    const Cluster *pBestMatchedCluster(nullptr);
    this->GetBestMatchedCluster(pCluster1, pCluster2, commonMuonPfoList, matchedClusterList, pBestMatchedCluster, reducedChiSquared);

    /*
    if (!pBestMatchedCluster)
    {
        return STATUS_CODE_NOT_FOUND;
    }
    */

    //ATTN: Ignore if other clusters matches have more hits
    if (pBestMatchedCluster && (pBestMatchedCluster->IsAvailable()))
    {
        unsigned int hitSum12(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits());
        unsigned int hitSum13(pCluster1->GetNCaloHits() + pBestMatchedCluster->GetNCaloHits());
        unsigned int hitSum23(pCluster2->GetNCaloHits() + pBestMatchedCluster->GetNCaloHits());        

        if (hitSum23 > hitSum12)
            return STATUS_CODE_NOT_FOUND;

        if (hitSum13 > hitSum12)
            return STATUS_CODE_NOT_FOUND;
    }

    TwoViewXOverlap xOverlapObject(xMin1, xMax1, xMin2, xMax2);

    overlapResult = TrackTwoViewTopologyOverlapResult(xOverlapObject, commonMuonPfoList, pBestMatchedCluster, matchedClusterList, reducedChiSquared);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::AreClustersCompatible(const Cluster *const pCluster1, const Cluster *const pCluster2, PfoList &commonMuonPfoList) const
{
    ClusterList consideredClusters1, consideredClusters2;
    PfoList nearbyMuonPfos1, nearbyMuonPfos2;
    this->GetNearbyMuonPfos(pCluster1, consideredClusters1, nearbyMuonPfos1);
    this->GetNearbyMuonPfos(pCluster2, consideredClusters2, nearbyMuonPfos2);

    for (const ParticleFlowObject *const pNearbyMuon1 : nearbyMuonPfos1)
    {
        for (const ParticleFlowObject *const pNearbyMuon2 : nearbyMuonPfos2)
        {
            if (pNearbyMuon1 == pNearbyMuon2)
            {
                commonMuonPfoList.push_back(pNearbyMuon1);
            }
        }
    }
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

float TwoViewDeltaRayMatchingAlgorithm::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        const float separationSquared((testPosition - referencePoint).GetMagnitudeSquared());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return std::sqrt(shortestDistanceSquared);
}


//------------------------------------------------------------------------------------------------------------------------------------------
    
void TwoViewDeltaRayMatchingAlgorithm::CollectHits(const Cluster *const pCluster1, const Cluster *const pCluster2, const CartesianPointVector &projectedPositions,
    ClusterList &matchedClusters) const
{
    const ClusterList *pInputClusterList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));
    
    HitOwnershipMap hitOwnershipMap;
    for (const Cluster *const pCluster : *pInputClusterList)
    {        
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            float separation(this->GetClosestDistance(pCaloHit, projectedPositions));
            
            if (separation < 2.f)
                hitOwnershipMap[pCluster].push_back(pCaloHit);
        }
    }

    for (auto &entry : hitOwnershipMap)
    {
        const Cluster *const pMatchedCluster(entry.first);

        XOverlap xThreeViewOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
        float chiSquaredSum(0.f);
        unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

        CaloHitList caloHitList1, caloHitList2, caloHitList3;
        pCluster1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
        pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
        pMatchedCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);
    
        StatusCode status(this->PerformMatching(caloHitList1, caloHitList2, caloHitList3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xThreeViewOverlapObject));
                        
        if (status == STATUS_CODE_NOT_FOUND)
            continue;

        if (status != STATUS_CODE_SUCCESS)
            throw StatusCodeException(status);

        const float reducedChiSquared(chiSquaredSum / nSamplingPoints);

        if (reducedChiSquared < 1.f)
            matchedClusters.push_back(entry.first);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GetBestMatchedCluster(const Cluster *const pCluster1, const Cluster *const pCluster2, const PfoList &commonMuonPfoList,
    const ClusterList &matchedClusters, const Cluster *&pBestMatchedCluster, float &reducedChiSquared) const
{
    const ClusterList *pInputClusterList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));
    
    ClusterList muonClusterList;
    
    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
        LArPfoHelper::GetClusters(pMuonPfo, LArClusterHelper::GetClusterHitType(pInputClusterList->front()), muonClusterList);
    
    unsigned int highestNHits(0);
    for (const Cluster *const pMatchedCluster : matchedClusters)
    {
        //std::cout << "pMatchedCluster: " << pMatchedCluster << std::endl;
        if (!pMatchedCluster->IsAvailable())
        {
            if (std::find(muonClusterList.begin(), muonClusterList.end(), pMatchedCluster) == muonClusterList.end())
                continue;
        }
                
        if (pMatchedCluster->GetNCaloHits() > highestNHits)
        {
            highestNHits = pMatchedCluster->GetNCaloHits();
            pBestMatchedCluster = pMatchedCluster;
        }
    }

    //std::cout << "left loop" << std::endl;
    
    if (!pBestMatchedCluster)
    {
        //std::cout << "left funciton: null ptr" << std::endl;
        return;
    }
    
    XOverlap xThreeViewOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    float chiSquaredSum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    CaloHitList caloHitList1, caloHitList2, caloHitList3;
    pCluster1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
    pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
    pBestMatchedCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);
    
    StatusCode status(this->PerformMatching(caloHitList1, caloHitList2, caloHitList3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xThreeViewOverlapObject));

    if (status == STATUS_CODE_NOT_FOUND)
    {
        std::cout << "ISOBEL: THIS SHOULD NEVER HAPPEN" << std::endl;
        throw StatusCodeException(status);
    }

    if (status != STATUS_CODE_SUCCESS)
    {
        std::cout << "here" << std::endl;
        throw StatusCodeException(status);
    }

    reducedChiSquared = (chiSquaredSum / nSamplingPoints);

    //std::cout << "left function: at end" << std::endl;
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

StatusCode TwoViewDeltaRayMatchingAlgorithm::PerformMatching(const CaloHitList &clusterU, const CaloHitList &clusterV, const CaloHitList &clusterW,
    float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &xOverlapObject) const
{    
    float xMinU(-std::numeric_limits<float>::max()), xMaxU(+std::numeric_limits<float>::max());
    float xMinV(-std::numeric_limits<float>::max()), xMaxV(+std::numeric_limits<float>::max());
    float xMinW(-std::numeric_limits<float>::max()), xMaxW(+std::numeric_limits<float>::max());

    this->GetClusterSpanX(clusterU, xMinU, xMaxU);
    this->GetClusterSpanX(clusterV, xMinV, xMaxV);
    this->GetClusterSpanX(clusterW, xMinW, xMaxW);

    // Need to remove the xPitch from calculations to be consistent with view xSpan calculated in the xOverlapObject
    const float xMinCentre(std::max(xMinU, std::max(xMinV, xMinW)));
    const float xMaxCentre(std::min(xMaxU, std::min(xMaxV, xMaxW)));
    const float xCentreOverlap(xMaxCentre - xMinCentre);

    if (xCentreOverlap < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;
    
    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMinU, std::max(xMinV, xMinW)) - xPitch);
    const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)) + xPitch);
    const float xOverlap(xMax - xMin);

    const HitType hitTypeU(clusterU.front()->GetHitType());
    const HitType hitTypeV(clusterV.front()->GetHitType());
    const HitType hitTypeW(clusterW.front()->GetHitType());

    if (hitTypeU == hitTypeV ||  hitTypeU == hitTypeW || hitTypeV == hitTypeW)
       return STATUS_CODE_FAILURE;

    const unsigned int nPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    chiSquaredSum = 0.f; nSamplingPoints = 0; nMatchedSamplingPoints = 0;
    
    for (unsigned int n = 0; n < nPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMinU(0.f), zMinV(0.f), zMinW(0.f), zMaxU(0.f), zMaxV(0.f), zMaxW(0.f);
            this->GetClusterSpanZ(clusterU, xmin, xmax, zMinU, zMaxU);
            this->GetClusterSpanZ(clusterV, xmin, xmax, zMinV, zMaxV);
            this->GetClusterSpanZ(clusterW, xmin, xmax, zMinW, zMaxW);

            const float zU(0.5f * (zMinU + zMaxU));
            const float zV(0.5f * (zMinV + zMaxV));
            const float zW(0.5f * (zMinW + zMaxW));

            const float dzU(zMaxU - zMinU);
            const float dzV(zMaxV - zMinV);
            const float dzW(zMaxW - zMinW);
            const float dzPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

            const float zprojU(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitTypeV, hitTypeW, zV, zW));
            const float zprojV(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitTypeW, hitTypeU, zW, zU));
            const float zprojW(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitTypeU, hitTypeV, zU, zV));

            ++nSamplingPoints;

            const float deltaSquared(((zU - zprojU) * (zU - zprojU) + (zV - zprojV) * (zV - zprojV) + (zW - zprojW) * (zW - zprojW)) / 3.f);
            const float sigmaSquared(dzU * dzU + dzV * dzV + dzW * dzW + dzPitch * dzPitch);
            const float pseudoChi2(deltaSquared / sigmaSquared);

            chiSquaredSum += pseudoChi2;
            
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

    // Apply tensor threshold cuts
    if (nSamplingPoints == 0)
        return STATUS_CODE_NOT_FOUND;

    const float matchedFraction(static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints));

    if ((matchedFraction < m_minMatchedFraction) || (nMatchedSamplingPoints < m_minMatchedPoints))
        return STATUS_CODE_NOT_FOUND;

    xOverlapObject = XOverlap(xMinU, xMaxU, xMinV, xMaxV, xMinW, xMaxW, xCentreOverlap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GetClusterSpanX(const CaloHitList &caloHitList, float &xMin, float &xMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    
    for(const CaloHit *const pCaloHit : caloHitList)
    {
        const float xCoordinate(pCaloHit->GetPositionVector().GetX());

        if (xCoordinate < xMin)
            xMin = xCoordinate;

        if (xCoordinate > xMax)
            xMax = xCoordinate;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::GetClusterSpanZ(const CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const
{
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();

    bool found(false);
    for(const CaloHit *const pCaloHit : caloHitList)
    {
        const float xCoordinate(pCaloHit->GetPositionVector().GetX());
        const float zCoordinate(pCaloHit->GetPositionVector().GetZ());
        
        if ((xCoordinate < xMin) || (xCoordinate > xMax))
            continue;

        found = true;

        if (zCoordinate < zMin)
            zMin = zCoordinate;

        if (zCoordinate > zMax)
            zMax = zCoordinate;
    }

    if (!found)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

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
