/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatchingAlgorithm.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

template<typename T>
NViewDeltaRayMatchingAlgorithm<T>::NViewDeltaRayMatchingAlgorithm() :
    m_searchRegion1D(3.f),    
    m_pseudoChi2Cut(1.5f), 
    m_xOverlapWindow(1.f),
    m_minMatchedFraction(0.5),
    m_minMatchedPoints(3),
    m_minProjectedPositions(3),
    m_maxCosmicRayHitFraction(0.05f),
    m_strayClusterSeparation(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (const Cluster *const pCluster : *pInputClusterList)
    {
        if ((pCluster->IsAvailable()) && (this->DoesClusterPassTensorThreshold(pCluster)))
            selectedClusterList.push_back(pCluster);
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::PrepareInputClusters(ClusterList &preparedClusterList)
{
    if (preparedClusterList.empty())
        return;

    const HitType &hitType(LArClusterHelper::GetClusterHitType(preparedClusterList.front()));

    this->FillHitToClusterMap(hitType);
    this->FillClusterProximityMap(hitType);
    this->FillClusterToPfoMap(hitType);
    this->FillStrayClusterList(hitType);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::FillHitToClusterMap(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));

    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterMap(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::AddToClusterMap(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
        hitToClusterMap[pCaloHit] = pCluster;
}     
        
//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::FillClusterProximityMap(const HitType &hitType)
{
    this->BuildKDTree(hitType);
    
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    
    for (const Cluster *const pCluster : inputClusterList)
        this->AddToClusterProximityMap(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::BuildKDTree(const HitType &hitType)
{
    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);    

    CaloHitList allCaloHits;

    for (auto &entry : hitToClusterMap)
        allCaloHits.push_back(entry.first);
    
    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));

    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::AddToClusterProximityMap(const Cluster *const pCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
  
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {             
        HitKDNode2DList found;
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));

        kdTree.search(searchRegionHits, found);
            
        for (const auto &hit : found)
        {
            const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));
            
            if (pNearbyCluster == pCluster)
                continue;
                
            ClusterList  &nearbyClusterList(clusterProximityMap[pCluster]);

            if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
                nearbyClusterList.push_back(pNearbyCluster);
            
            ClusterList &invertedNearbyClusterList(clusterProximityMap[pNearbyCluster]);
            
            if (std::find(invertedNearbyClusterList.begin(), invertedNearbyClusterList.end(), pCluster) == invertedNearbyClusterList.end())
                invertedNearbyClusterList.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------        

template<typename T>    
void NViewDeltaRayMatchingAlgorithm<T>::FillClusterToPfoMap(const HitType &hitType)
{    
    const PfoList *pMuonPfoList(nullptr);
    
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_muonPfoListName, pMuonPfoList));

    if ((!pMuonPfoList) || pMuonPfoList->empty())
        return;

    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);      

    for (const ParticleFlowObject *const pPfo : *pMuonPfoList)
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetClusters(pPfo, hitType, pfoClusters);

        for (const Cluster *const pCluster : pfoClusters)
            clusterToPfoMap[pCluster] = pPfo;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::FillStrayClusterList(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);

    for (const Cluster *const pCluster : inputClusterList)
    {
        if ((pCluster->IsAvailable()) && (!this->DoesClusterPassTensorThreshold(pCluster)))
            strayClusterList.push_back(pCluster);
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::GetMuonCluster(const PfoList &commonMuonPfoList, const HitType &hitType, const Cluster *&pMuonCluster) const
{
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

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::GetNearbyMuonPfos(const Cluster *const pCluster, ClusterList &consideredClusters, PfoList &nearbyMuonPfos) const 
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
    const ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);

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

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::PerformThreeViewMatching(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const pCluster3,
    float &reducedChiSquared) const
{
    float chiSquaredSum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    XOverlap xThreeViewOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

    if (this->PerformThreeViewMatching(pCluster1, pCluster2, pCluster3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xThreeViewOverlapObject) == STATUS_CODE_NOT_FOUND)
        return STATUS_CODE_NOT_FOUND;

    reducedChiSquared = chiSquaredSum / nSamplingPoints;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::PerformThreeViewMatching(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &xOverlapObject) const
{    
    float xMinU(-std::numeric_limits<float>::max()), xMaxU(+std::numeric_limits<float>::max());
    float xMinV(-std::numeric_limits<float>::max()), xMaxV(+std::numeric_limits<float>::max());
    float xMinW(-std::numeric_limits<float>::max()), xMaxW(+std::numeric_limits<float>::max());

    pClusterU->GetClusterSpanX(xMinU, xMaxU);
    pClusterV->GetClusterSpanX(xMinV, xMaxV);
    pClusterW->GetClusterSpanX(xMinW, xMaxW);

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

    const HitType hitTypeU(LArClusterHelper::GetClusterHitType(pClusterU));
    const HitType hitTypeV(LArClusterHelper::GetClusterHitType(pClusterV));
    const HitType hitTypeW(LArClusterHelper::GetClusterHitType(pClusterW));

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
            pClusterU->GetClusterSpanZ(xmin, xmax, zMinU, zMaxU);
            pClusterV->GetClusterSpanZ(xmin, xmax, zMinV, zMaxV);
            pClusterW->GetClusterSpanZ(xmin, xmax, zMinW, zMaxW);

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

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::PerformThreeViewMatching(const CaloHitList &pCluster1, const CaloHitList &pCluster2, const CaloHitList &pCluster3,
    float &reducedChiSquared) const
{
    float chiSquaredSum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    XOverlap xThreeViewOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

    if (this->PerformThreeViewMatching(pCluster1, pCluster2, pCluster3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xThreeViewOverlapObject) == STATUS_CODE_NOT_FOUND)
        return STATUS_CODE_NOT_FOUND;

    reducedChiSquared = chiSquaredSum / nSamplingPoints;

    return STATUS_CODE_SUCCESS;
}
   
//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::PerformThreeViewMatching(const CaloHitList &clusterU, const CaloHitList &clusterV, const CaloHitList &clusterW,
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

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::GetClusterSpanX(const CaloHitList &caloHitList, float &xMin, float &xMax) const
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

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::GetClusterSpanZ(const CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const
{
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();

    bool found(false);
    for (const CaloHit *const pCaloHit : caloHitList)
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

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::ProjectMuonPositions(const HitType &thirdViewHitType, const ParticleFlowObject *const pParentMuon,
    CartesianPointVector &projectedPositions) const
{
    ClusterList muonClusterList1, muonClusterList2;    

    for (const HitType &hitType1 : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        if (hitType1 == thirdViewHitType)
            continue;
        
        for (const HitType &hitType2 : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if ((hitType2 == thirdViewHitType) || (hitType1 == hitType2))
                continue;

            LArPfoHelper::GetClusters(pParentMuon, hitType1, muonClusterList1);
            LArPfoHelper::GetClusters(pParentMuon, hitType2, muonClusterList2);
            
            if ((muonClusterList1.size() != 1) || (muonClusterList2.size() != 1))
                return STATUS_CODE_NOT_FOUND;
            
            break;
        }

        break;
    }

    return (this->GetProjectedPositions(muonClusterList1.front(), muonClusterList2.front(), projectedPositions));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::GetProjectedPositions(const Cluster *const pCluster1, const Cluster *const pCluster2,
    CartesianPointVector &projectedPositions) const
{
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1); 
    pCluster2->GetClusterSpanX(xMin2, xMax2); 

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, xMin2) - xPitch);
    const float xMax(std::min(xMax1, xMax2) + xPitch);
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
         return STATUS_CODE_NOT_FOUND;
    
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    if (hitType1 == hitType2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

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

    // Reject if projection is not good
    if (projectedPositions.size() < m_minProjectedPositions)
        return STATUS_CODE_NOT_FOUND;
    
    return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::CollectHitsFromMuon(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const Cluster *const pThirdViewCluster, const ParticleFlowObject *const pParentMuon, const float minDistanceFromMuon, const float maxDistanceToCollected,
    CaloHitList &collectedHits) const
{
    HitType thirdViewHitType(TPC_VIEW_U);
    CartesianPointVector deltaRayProjectedPositions;
    
    if (!pThirdViewCluster)
    {
        if (this->GetProjectedPositions(pCluster1, pCluster2, deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
            return STATUS_CODE_NOT_FOUND;

        for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if ((LArClusterHelper::GetClusterHitType(pCluster1) != hitType) && (LArClusterHelper::GetClusterHitType(pCluster2) != hitType))
            {
                thirdViewHitType = hitType;
                break;
            }
        }
    }
    else
    {
        CaloHitList deltaRayCaloHitList;
        pThirdViewCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayCaloHitList);

        for (const CaloHit *const pCaloHit : deltaRayCaloHitList)
            deltaRayProjectedPositions.push_back(pCaloHit->GetPositionVector());

        thirdViewHitType = LArClusterHelper::GetClusterHitType(pThirdViewCluster);
    }

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pParentMuon, thirdViewHitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    const Cluster *const pMuonCluster(muonClusterList.front());

    // To avoid fluctuatuions, parameterise the muon track
    CartesianVector positionOnMuon(0.f, 0.f, 0.f), muonDirection(0.f, 0.f, 0.f);
    if (this->ParameteriseMuon(pParentMuon, deltaRayProjectedPositions, thirdViewHitType, positionOnMuon, muonDirection) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;
    
    this->CollectHitsFromMuon(positionOnMuon, muonDirection, pMuonCluster, deltaRayProjectedPositions, minDistanceFromMuon, maxDistanceToCollected, collectedHits);

    if (collectedHits.empty())
        return STATUS_CODE_NOT_FOUND;

    // Catch if delta ray has travelled along muon
    if ((static_cast<float>(collectedHits.size()) / pMuonCluster->GetNCaloHits()) > m_maxCosmicRayHitFraction)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::ParameteriseMuon(const ParticleFlowObject *const pParentMuon, const Cluster *const pDeltaRayCluster,
    CartesianVector &positionOnMuon, CartesianVector &muonDirection) const
{
    CaloHitList deltaRayHitList;
    pDeltaRayCluster->GetOrderedCaloHitList().FillCaloHitList(deltaRayHitList);

    CartesianPointVector deltaRayProjectedPositions;

    for (const CaloHit *const pCaloHit : deltaRayHitList)
        deltaRayProjectedPositions.push_back(pCaloHit->GetPositionVector());

    return this->ParameteriseMuon(pParentMuon, deltaRayProjectedPositions, LArClusterHelper::GetClusterHitType(pDeltaRayCluster), positionOnMuon, muonDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::ParameteriseMuon(const ParticleFlowObject *const pParentMuon, const CartesianPointVector &deltaRayProjectedPositions,
    const HitType &hitType, CartesianVector &positionOnMuon, CartesianVector &muonDirection) const
{
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pParentMuon, hitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return STATUS_CODE_NOT_FOUND;
    
    CartesianPointVector muonProjectedPositions;
    if (this->ProjectMuonPositions(hitType, pParentMuon, muonProjectedPositions) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    const Cluster *const pMuonCluster(muonClusterList.front());
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 40, slidingFitPitch);

    CartesianVector deltaRayVertex(0.f, 0.f, 0.f), muonVertex(0.f, 0.f, 0.f);
    LArMuonLeadingHelper::GetClosestPositions(deltaRayProjectedPositions, pMuonCluster, deltaRayVertex, muonVertex);

    positionOnMuon = LArMuonLeadingHelper::GetClosestPosition(muonVertex, muonProjectedPositions, pMuonCluster);

    if (positionOnMuon.GetMagnitude() < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(positionOnMuon, rL, rT);
    slidingFitResult.GetGlobalFitDirection(rL, muonDirection);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::CollectHitsFromMuon(const CartesianVector &positionOnMuon, const CartesianVector &muonDirection,
    const Cluster *const pMuonCluster, const CartesianPointVector &deltaRayProjectedPositions, const float &minDistanceFromMuon, const float maxDistanceToCollected,
     CaloHitList &collectedHits) const
{
    CaloHitList cosmicRayHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(cosmicRayHitList);
    
    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

	    for (const CaloHit *const pCaloHit : cosmicRayHitList)
        {
            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const float distanceToCollectedHits(std::min(LArMuonLeadingHelper::GetClosestDistance(pCaloHit, deltaRayProjectedPositions),
                LArMuonLeadingHelper::GetClosestDistance(pCaloHit, collectedHits)));
            const float distanceToMuonHits(muonDirection.GetCrossProduct(pCaloHit->GetPositionVector() - positionOnMuon).GetMagnitude());

            if ((std::fabs(distanceToMuonHits - distanceToCollectedHits) < std::numeric_limits<float>::epsilon()) || (distanceToMuonHits < minDistanceFromMuon) ||
                (distanceToCollectedHits > distanceToMuonHits) || (distanceToCollectedHits > maxDistanceToCollected))
            {
                continue;
            }
                
            collectedHits.push_back(pCaloHit);
            hitsAdded = true;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::SplitMuonCluster(const std::string &clusterListName, const Cluster *const pMuonCluster,
    const CaloHitList &collectedHits, const Cluster *&pDeltaRayCluster) const
{
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

    CaloHitList muonCaloHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);
    
    for (const CaloHit *const pCaloHit : muonCaloHitList)
    {
        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pMuonCluster, pCaloHit));

            if (!pDeltaRayCluster)
            {
                const ClusterList *pTemporaryList(nullptr);
                std::string temporaryListName, currentListName;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pTemporaryList, temporaryListName));

                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pDeltaRayCluster));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pDeltaRayCluster, pCaloHit));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool NViewDeltaRayMatchingAlgorithm<T>::CreatePfos(ProtoParticleVector &protoParticleVector)
{
    std::sort(protoParticleVector.begin(), protoParticleVector.end(), [] (const ProtoParticle &a, const ProtoParticle &b) -> bool
    {
        unsigned int aHitTotal(0);
        for (const Cluster *const pCluster : a.m_clusterList)
            aHitTotal += pCluster->GetNCaloHits();

        unsigned int bHitTotal(0);
        for (const Cluster *const pCluster : b.m_clusterList)
            bHitTotal += pCluster->GetNCaloHits();

        return (aHitTotal > bHitTotal);
    });

    for (ProtoParticle protoParticle : protoParticleVector)
    {
        float longestSpan(-std::numeric_limits<float>::max()), spanMinX(0.f), spanMaxX(0.f);
    
        for (const Cluster *const pCluster : protoParticle.m_clusterList)
        {
            float minX(0.f), maxX(0.f);
            pCluster->GetClusterSpanX(minX, maxX);

            const float span(maxX - minX);
            
            if (span > longestSpan)
            {
                longestSpan = span; spanMinX = minX; spanMaxX = maxX;
            }
        }

	    for (const Cluster *const pCluster : protoParticle.m_clusterList)
	    {            
            ClusterList collectedClusters;
            this->CollectStrayClusters(pCluster, spanMinX, spanMaxX, collectedClusters);

	         if (!collectedClusters.empty())
                 this->AddInStrayClusters(pCluster, collectedClusters);
	    }
	}

    return (this->CreateThreeDParticles(protoParticleVector));
}      

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::CollectStrayClusters(const Cluster *const pClusterToEnlarge, const float rangeMinX, const float rangeMaxX, ClusterList &collectedClusters)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pClusterToEnlarge));
    const ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);
    const ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pClusterToEnlarge));

    if (clusterProximityIter == clusterProximityMap.end())
        return;

    for (const Cluster *const pNearbyCluster : clusterProximityIter->second)
    {
        if (std::find(strayClusterList.begin(), strayClusterList.end(), pNearbyCluster) == strayClusterList.end())
            continue;

        float xMin(-std::numeric_limits<float>::max()), xMax(+std::numeric_limits<float>::max());
        pNearbyCluster->GetClusterSpanX(xMin, xMax);

        if (!(((xMin > rangeMinX) && (xMin < rangeMaxX)) || ((xMax > rangeMinX) && (xMax < rangeMaxX))))
            continue;
        
        if (LArClusterHelper::GetClosestDistance(pClusterToEnlarge, pNearbyCluster) > m_strayClusterSeparation)
            continue;

        if (std::find(collectedClusters.begin(), collectedClusters.end(), pNearbyCluster) == collectedClusters.end())
            collectedClusters.push_back(pNearbyCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::AddInStrayClusters(const Cluster *const pClusterToEnlarge, const ClusterList &collectedClusters)
{
    this->UpdateUponDeletion(pClusterToEnlarge);

    for(const Cluster *const pCollectedCluster : collectedClusters)
    {
        this->UpdateUponDeletion(pCollectedCluster);
        
        std::string clusterListName(this->GetClusterListName(LArClusterHelper::GetClusterHitType(pClusterToEnlarge)));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pCollectedCluster, clusterListName, clusterListName));
    }

    this->UpdateContainers({pClusterToEnlarge}, {nullptr});
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pDeletedCluster));
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);    
    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
    ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);

    CaloHitList caloHitList;
    pDeletedCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const HitToClusterMap::const_iterator iter(hitToClusterMap.find(pCaloHit));

        if (iter == hitToClusterMap.end())
            throw StatusCodeException(STATUS_CODE_FAILURE);
        
        hitToClusterMap.erase(iter);
    }

    const ClusterList::const_iterator strayClusterIter(std::find(strayClusterList.begin(), strayClusterList.end(), pDeletedCluster));

    if (strayClusterIter != strayClusterList.end())
        strayClusterList.erase(strayClusterIter);    

    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pDeletedCluster));

    if (clusterProximityIter != clusterProximityMap.end())
    {
        const ClusterList &nearbyClusterList(clusterProximityIter->second);

        for (const Cluster *const pNearbyCluster : nearbyClusterList)
        {
            const ClusterProximityMap::iterator iter(clusterProximityMap.find(pNearbyCluster));

            if (iter == clusterProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);

            ClusterList::iterator invertedIter(std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pDeletedCluster));
            invertedCloseClusters.erase(invertedIter);
        }
        
        clusterProximityMap.erase(clusterProximityIter);
    }

    const ClusterToPfoMap::const_iterator clusterToPfoIter(clusterToPfoMap.find(pDeletedCluster));

    if (clusterToPfoIter != clusterToPfoMap.end())
    {
        clusterToPfoMap.erase(clusterToPfoIter);
    }
    else
    {
        NViewMatchingAlgorithm<T>::UpdateUponDeletion(pDeletedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::UpdateForNewClusters(const ClusterVector &newClusterVector, const PfoVector &pfoVector)
{
    this->UpdateContainers(newClusterVector, pfoVector);
    
    for (unsigned int i = 0; i < newClusterVector.size(); i++)
    {
        const Cluster *const pNewCluster(newClusterVector.at(i));
        const ParticleFlowObject *const pMuonPfo(pfoVector.at(i));

        // ATTN: Only add delta ray clusters into the tensor
        if (!pMuonPfo)
            NViewMatchingAlgorithm<T>::UpdateForNewCluster(pNewCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::UpdateContainers(const ClusterVector &newClusterVector, const PfoVector &pfoVector)
{
    for (const Cluster *const pNewCluster : newClusterVector)
        this->AddToClusterMap(pNewCluster);
    
    for (unsigned int i = 0; i < newClusterVector.size(); i++)
    {
        const Cluster *const pNewCluster(newClusterVector.at(i));
        const ParticleFlowObject *const pMuonPfo(pfoVector.at(i));
        
        this->AddToClusterProximityMap(pNewCluster);
        
        if (pMuonPfo)
        {
            const HitType &hitType(LArClusterHelper::GetClusterHitType(pNewCluster));
            ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
            clusterToPfoMap[pNewCluster] = pMuonPfo;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void NViewDeltaRayMatchingAlgorithm<T>::TidyUp()
{
    m_hitToClusterMapU.clear();
    m_hitToClusterMapV.clear();
    m_hitToClusterMapW.clear();    

    m_kdTreeU.clear();
    m_kdTreeV.clear();
    m_kdTreeW.clear();    

    m_clusterProximityMapU.clear();
    m_clusterProximityMapV.clear();
    m_clusterProximityMapW.clear();
    
    m_clusterToPfoMapU.clear();
    m_clusterToPfoMapV.clear();
    m_clusterToPfoMapW.clear();

    m_strayClusterListU.clear();
    m_strayClusterListV.clear();
    m_strayClusterListW.clear();
    
    return NViewMatchingAlgorithm<T>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode NViewDeltaRayMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName)); 
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SearchRegion1D", m_searchRegion1D));        

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverlapWindow", m_xOverlapWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedPoints", m_minMatchedPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinProjectedPositions", m_minProjectedPositions));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCosmicRayHitFraction", m_maxCosmicRayHitFraction));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StrayClusterSeparation", m_strayClusterSeparation));    
    
    return NViewMatchingAlgorithm<T>::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class NViewDeltaRayMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult>>;
template class NViewDeltaRayMatchingAlgorithm<TwoViewMatchingControl<TwoViewDeltaRayOverlapResult>>;

} // namespace lar_content

