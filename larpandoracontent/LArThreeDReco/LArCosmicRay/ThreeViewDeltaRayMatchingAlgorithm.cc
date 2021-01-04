/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

ThreeViewDeltaRayMatchingAlgorithm::ThreeViewDeltaRayMatchingAlgorithm()  :
    m_isStrayListUInitialised(false),
    m_isStrayListVInitialised(false),
    m_isStrayListWInitialised(false),
    m_strayClusterListU(ClusterList()),
    m_strayClusterListV(ClusterList()),
    m_strayClusterListW(ClusterList()),
    m_nMaxTensorToolRepeats(10),
    m_minClusterCaloHits(5),//5),
    m_searchRegion1D(3.f),    
    m_pseudoChi2Cut(1.5f), //normally three
    m_xOverlapWindow(1.f),
    m_minMatchedFraction(0.5),
    m_minMatchedPoints(3)//3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (const Cluster *const pCluster : *pInputClusterList)
    {
        if ((pCluster->IsAvailable()) && (this->DoesClusterPassTesorThreshold(pCluster)))
            selectedClusterList.push_back(pCluster);
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewDeltaRayMatchingAlgorithm::DoesClusterPassTesorThreshold(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
        return false;

    return true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::PrepareInputClusters(ClusterList &preparedClusterList)
{
    if (preparedClusterList.empty())
        return;

    const HitType &hitType(LArClusterHelper::GetClusterHitType(preparedClusterList.front()));
    
    this->FillHitToClusterMap(hitType);
    this->FillClusterProximityMap(hitType);
    this->FillClusterToPfoMap(hitType);
    this->InitialiseStrayClusterList(hitType);
}
        
//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    DeltaRayOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult));

    if (overlapResult.IsInitialized())
    {
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    DeltaRayOverlapResult &overlapResult) const
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
    
    PfoList commonMuonPfoList;
    this->AreClustersCompatible(pClusterU, pClusterV, pClusterW, commonMuonPfoList);
    
    if (commonMuonPfoList.empty())
        return STATUS_CODE_NOT_FOUND;
    
    // what is the m_xOverlapWindow? - seems to be like a hit width? (an uncertainty in x)
    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMinU, std::max(xMinV, xMinW)) - xPitch);
    const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)) + xPitch);
    const float xOverlap(xMax - xMin);

    const HitType hitTypeU(LArClusterHelper::GetClusterHitType(pClusterU));
    const HitType hitTypeV(LArClusterHelper::GetClusterHitType(pClusterV));
    const HitType hitTypeW(LArClusterHelper::GetClusterHitType(pClusterW));

    if (hitTypeU != TPC_VIEW_U ||  hitTypeV != TPC_VIEW_V || hitTypeW != TPC_VIEW_W)
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
    
    // Apply tensor threshold cuts
    if (nSamplingPoints == 0)
        return STATUS_CODE_NOT_FOUND;

    const float matchedFraction(static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints));

    if ((matchedFraction < m_minMatchedFraction) || (nMatchedSamplingPoints < m_minMatchedPoints))
        return STATUS_CODE_NOT_FOUND;

    const XOverlap xOverlapObject(xMinU, xMaxU, xMinV, xMaxV, xMinW, xMaxW, xCentreOverlap);

    overlapResult = DeltaRayOverlapResult(nMatchedSamplingPoints, nSamplingPoints, pseudoChi2Sum, xOverlapObject, commonMuonPfoList);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::AreClustersCompatible(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    PfoList &commonMuonPfoList) const
{
    ClusterList consideredClustersU, consideredClustersV, consideredClustersW;
    PfoList nearbyMuonPfosU, nearbyMuonPfosV, nearbyMuonPfosW;
    this->GetNearbyMuonPfos(pClusterU, consideredClustersU, nearbyMuonPfosU);
    this->GetNearbyMuonPfos(pClusterV, consideredClustersV, nearbyMuonPfosV);
    this->GetNearbyMuonPfos(pClusterW, consideredClustersW, nearbyMuonPfosW);

    for (const ParticleFlowObject *const pNearbyMuonU : nearbyMuonPfosU)
    {
        for (const ParticleFlowObject *const pNearbyMuonV : nearbyMuonPfosV)
        {
            for (const ParticleFlowObject *const pNearbyMuonW : nearbyMuonPfosW)
            {
                if ((pNearbyMuonU == pNearbyMuonV) && (pNearbyMuonV == pNearbyMuonW))
                {
                    commonMuonPfoList.push_back(pNearbyMuonU);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::GetNearbyMuonPfos(const Cluster *const pCluster, ClusterList &consideredClusters, PfoList &nearbyMuonPfos) const 
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

void ThreeViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (TensorToolVector::const_iterator toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end(); )
    {
        if ((*toolIter)->Run(this, this->GetMatchingControl().GetOverlapTensor()))
        {
            toolIter = m_algorithmToolVector.begin();

	    //std::cout << "tool counter: " << repeatCounter << std::endl;

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++toolIter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::TidyUp()
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

    this->ClearStrayClusterLists();
    
    return BaseAlgorithm::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ClusterRebuilding", m_reclusteringAlgorithmName));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName));
    
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "DeltaRayTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DeltaRayTensorTool *const pDeltaRayTensorTool(dynamic_cast<DeltaRayTensorTool*>(*iter));

        if (!pDeltaRayTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayTensorTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

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
    
    return BaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::FillHitToClusterMap(const HitType &hitType)
{
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);

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

void ThreeViewDeltaRayMatchingAlgorithm::FillClusterProximityMap(const HitType &hitType)
{
    const HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);

    CaloHitList allCaloHits;
    for (auto &entry : hitToClusterMap)
        allCaloHits.push_back(entry.first);
    
    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));    
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);    
    
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

void ThreeViewDeltaRayMatchingAlgorithm::FillClusterToPfoMap(const HitType &hitType)
{
    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);        
    
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

void ThreeViewDeltaRayMatchingAlgorithm::UpdateForNewClusters(const ClusterVector &newClusterVector, const PfoVector &pfoVector)
{
    for (const Cluster *const pNewCluster : newClusterVector)
    {
        const HitType &hitType(LArClusterHelper::GetClusterHitType(pNewCluster));
	HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);

        CaloHitList caloHitList;
	pNewCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *const pCaloHit : caloHitList)
        { 
            hitToClusterMap[pCaloHit] = pNewCluster;
	}
    }
 
   for (unsigned int i = 0; i < newClusterVector.size(); i++)
    {
        const Cluster *const pNewCluster(newClusterVector.at(i));

        const HitType &hitType(LArClusterHelper::GetClusterHitType(pNewCluster));

	HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
	HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);
	ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
	ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);
  
        CaloHitList caloHitList;
	pNewCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *const pCaloHit : caloHitList)
        {             
            HitKDNode2DList found;  
	    KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));

	    kdTree.search(searchRegionHits, found);
            
	    for (const auto &hit : found)
            {
	        if (std::find(caloHitList.begin(), caloHitList.end(), hit.data) != caloHitList.end())
		  continue;
                
		ClusterList  &nearbyClusterList(clusterProximityMap[pNewCluster]);
		const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));

		if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
		  nearbyClusterList.push_back(pNearbyCluster);
	    }
	}

	if (clusterProximityMap.find(pNewCluster) != clusterProximityMap.end())
        {
	    const ClusterList &nearbyClusterList(clusterProximityMap.at(pNewCluster));

	    for (const Cluster *const pNearbyCluster : nearbyClusterList)
            {
	        ClusterList &invertedCloseClusters(clusterProximityMap[pNearbyCluster]);
            
		if (std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pNewCluster) == invertedCloseClusters.end())
		    invertedCloseClusters.push_back(pNewCluster);
	    }
	}

	const ParticleFlowObject *const pMuonPfo(pfoVector.at(i));

        if (pMuonPfo)
        {
            clusterToPfoMap[pNewCluster] = pMuonPfo;
	}
    }

    for (unsigned int i = 0; i < newClusterVector.size(); i++)
    {
        const Cluster *const pNewCluster(newClusterVector.at(i));
	const ParticleFlowObject *const pMuonPfo(pfoVector.at(i));

	if (!pMuonPfo)
	  BaseAlgorithm::UpdateForNewCluster(pNewCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::UpdateForCreation(const Cluster *const pNewCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pNewCluster));

    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);    
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);    

    CaloHitList caloHitList;
    pNewCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
        hitToClusterMap[pCaloHit] = pNewCluster;
  
    for (const CaloHit *const pCaloHit : caloHitList)
    {             
        HitKDNode2DList found;
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));

        kdTree.search(searchRegionHits, found);
            
        for (const auto &hit : found)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), hit.data) != caloHitList.end())
                continue;
                
            ClusterList  &nearbyClusterList(clusterProximityMap[pNewCluster]);
            const Cluster *const pNearbyCluster(hitToClusterMap.at(hit.data));

            if (std::find(nearbyClusterList.begin(), nearbyClusterList.end(), pNearbyCluster) == nearbyClusterList.end())
                nearbyClusterList.push_back(pNearbyCluster);
        }
    }
    
    if (clusterProximityMap.find(pNewCluster) != clusterProximityMap.end())
    {
        const ClusterList &nearbyClusterList(clusterProximityMap.at(pNewCluster));

        for (const Cluster *const pNearbyCluster : nearbyClusterList)
        {
            ClusterList &invertedCloseClusters(clusterProximityMap[pNearbyCluster]);
            
            if (std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pNewCluster) == invertedCloseClusters.end())
                invertedCloseClusters.push_back(pNewCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pDeletedCluster));

    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);
    ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

    const OrderedCaloHitList &orderedCaloHitList(pDeletedCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            auto iter(hitToClusterMap.find(pCaloHit));
            hitToClusterMap.erase(iter);
        }
    }

    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pDeletedCluster));

    if (clusterProximityIter != clusterProximityMap.end())
    {
        ClusterVector closeClusterVector;

        for (const Cluster *const pCloseCluster : clusterProximityIter->second)
            closeClusterVector.push_back(pCloseCluster);

        std::sort(closeClusterVector.begin(), closeClusterVector.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pCloseCluster : closeClusterVector)
        {
            const ClusterProximityMap::iterator iter(clusterProximityMap.find(pCloseCluster));

            if (iter == clusterProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);
            invertedCloseClusters.remove(pDeletedCluster);
        }
        clusterProximityMap.erase(clusterProximityIter);
    }

    const ClusterToPfoMap::const_iterator clusterToPfoIter(clusterToPfoMap.find(pDeletedCluster));

    if (clusterToPfoIter != clusterToPfoMap.end())
    {
      //std::cout << "MUON" << std::endl;
      clusterToPfoMap.erase(clusterToPfoIter);
    }
    else
    {
      //std::cout << "DELTA RAY" << std::endl;
      BaseAlgorithm::UpdateUponDeletion(pDeletedCluster);
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::InitialiseStrayClusterList(const HitType &hitType)
{
    if (this->IsStrayClusterListInitialised(hitType))
        return;

    auto &theTensor(this->GetMatchingControl().GetOverlapTensor());
    
    ClusterSet clusterSetU, clusterSetV, clusterSetW;
    theTensor.GetClustersInTensor(clusterSetU, clusterSetV, clusterSetW);

    ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);
    
    const ClusterList &inputClusterList(this->GetInputClusterList(hitType));

    for (const Cluster *const pCluster : inputClusterList)
    {
      if ((!this->DoesClusterPassTesorThreshold(pCluster)) && (pCluster->IsAvailable()) && (!theTensor.IsInTensor(pCluster, clusterSetU, clusterSetV, clusterSetW)))
      {
	if (pCluster->GetNCaloHits() > 5)
	{
	  std::cout << "TOO MANY HITS" << std::endl;
	}
            strayClusterList.push_back(pCluster);
      }
    }

    bool &isStrayClusterListInitialised((hitType == TPC_VIEW_U) ? m_isStrayListUInitialised : (hitType == TPC_VIEW_V) ? m_isStrayListVInitialised : m_isStrayListWInitialised);
    isStrayClusterListInitialised = true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewDeltaRayMatchingAlgorithm::IsStrayClusterListInitialised(const HitType &hitType) const
{
    if ((hitType != TPC_VIEW_U) && (hitType != TPC_VIEW_V) && (hitType != TPC_VIEW_W))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    return ((hitType == TPC_VIEW_U) ? m_isStrayListUInitialised : (hitType == TPC_VIEW_V) ? m_isStrayListVInitialised : m_isStrayListWInitialised);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::ClearStrayClusterLists()
{
    if (m_isStrayListUInitialised)
        m_strayClusterListU.clear();

    if (m_isStrayListVInitialised)
        m_strayClusterListV.clear();

    if (m_isStrayListWInitialised)
        m_strayClusterListW.clear();
    
    m_isStrayListUInitialised = false;
    m_isStrayListVInitialised = false;
    m_isStrayListWInitialised = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterList &ThreeViewDeltaRayMatchingAlgorithm::GetStrayClusterList(const HitType &hitType) const
{
    if ((hitType != TPC_VIEW_U) && (hitType != TPC_VIEW_V) && (hitType != TPC_VIEW_W))
        throw STATUS_CODE_NOT_ALLOWED;

    return ((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::RemoveFromStrayClusterList(const Cluster *const pClusterToRemove)
{
  const HitType &hitType(LArClusterHelper::GetClusterHitType(pClusterToRemove));

  if ((hitType != TPC_VIEW_U) && (hitType != TPC_VIEW_V) && (hitType != TPC_VIEW_W))
    throw STATUS_CODE_NOT_ALLOWED;

  ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);
  
  const ClusterList::const_iterator iter(std::find(strayClusterList.begin(), strayClusterList.end(), pClusterToRemove));

  if (iter == strayClusterList.end())
  {
    std::cout << "ISOBEL - NOT IN STRAY CLUSTER LIST" << std::endl;
    throw;
  }

  if (iter != strayClusterList.end())
    strayClusterList.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::CollectStrayHits(const Cluster *const pBadCluster, const float spanMinX, const float spanMaxX, ClusterList &collectedClusters)
{
    //std::cout << "AAAAAA" << std::endl;
    const HitType badHitType(LArClusterHelper::GetClusterHitType(pBadCluster));

    //std::cout << "BBBBBBB" << std::endl;
    
    this->InitialiseStrayClusterList(badHitType);

    //std::cout << "CCCCCCCC" << std::endl;
    //std::cout << "IS INITIALISED? " << this->IsStrayClusterListInitialised(badHitType) << std::endl;
    const ClusterList &strayClusterList(this->GetStrayClusterList(badHitType));
    //std::cout << "DDDDDDD" << std::endl;    
    //std::cout << "size: " << strayClusterList.size() << std::endl;

    for (const Cluster *const pCluster : strayClusterList)
    {
        float xMin(-std::numeric_limits<float>::max()), xMax(+std::numeric_limits<float>::max());
        //std::cout << "1111111" << std::endl;
	//std::cout << "cluster hit number: " << pCluster->GetNCaloHits() << std::endl;
        //std::cout << "JAM" << std::endl;
        pCluster->GetClusterSpanX(xMin, xMax);
        
        //std::cout << "22222222" << std::endl;        

        if (!(((xMin > spanMinX) && (xMin < spanMaxX)) || ((xMax > spanMinX) && (xMax < spanMaxX))))
            continue;

        //std::cout << "3333333" << std::endl;      
        
        if (LArClusterHelper::GetClosestDistance(pBadCluster, pCluster) > 2.f)
            continue;

        //std::cout << "4444444" << std::endl;    
	if (std::find(collectedClusters.begin(), collectedClusters.end(), pCluster) == collectedClusters.end())
	  collectedClusters.push_back(pCluster);
    }
    //std::cout << "EEEEEEEE" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void ThreeViewDeltaRayMatchingAlgorithm::AddInStrayClusters(const Cluster *const pClusterToEnlarge, const ClusterList &collectedClusters)
{
    this->UpdateUponDeletion(pClusterToEnlarge);

    for(const Cluster *const pCollectedCluster : collectedClusters)
    {
      //const ClusterList &strayClusterList1(this->GetStrayClusterList(LArClusterHelper::GetClusterHitType(pClusterToEnlarge)));
      //std::cout << "clusterList1 size: " << strayClusterList1.size() << std::endl;
        
        this->RemoveFromStrayClusterList(pCollectedCluster);

        //const ClusterList &strayClusterList2(this->GetStrayClusterList(LArClusterHelper::GetClusterHitType(pClusterToEnlarge)));
        //std::cout << "clusterList2 size: " << strayClusterList2.size() << std::endl;
        this->UpdateUponDeletion(pCollectedCluster);

        std::string clusterListName(this->GetClusterListName(LArClusterHelper::GetClusterHitType(pClusterToEnlarge)));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pCollectedCluster, clusterListName, clusterListName));
    }

    this->UpdateForCreation(pClusterToEnlarge);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeViewDeltaRayMatchingAlgorithm::CalculateChiSquared(const CaloHitList &clusterU, const CaloHitList &clusterV, const CaloHitList &clusterW) const
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
        return -1.f;//std::numeric_limits<float>::max();
    
    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMinU, std::max(xMinV, xMinW)) - xPitch);
    const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)) + xPitch);
    const float xOverlap(xMax - xMin);

    const HitType hitTypeU(clusterU.front()->GetHitType());
    const HitType hitTypeV(clusterV.front()->GetHitType());
    const HitType hitTypeW(clusterW.front()->GetHitType());

    //if (hitTypeU != TPC_VIEW_U ||  hitTypeV != TPC_VIEW_V || hitTypeW != TPC_VIEW_W)
    //return -2.f;//STATUS_CODE_FAILURE;

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

    if (nMatchedSamplingPoints < 3)
        return -3.f;//std::numeric_limits<float>::max();

    return (pseudoChi2Sum / static_cast<float>(nSamplingPoints));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::GetClusterSpanX(const CaloHitList &caloHitList, float &xMin, float &xMax) const
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

StatusCode ThreeViewDeltaRayMatchingAlgorithm::GetClusterSpanZ(const CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const
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

bool ThreeViewDeltaRayMatchingAlgorithm::CreatePfos(ProtoParticleVector &protoParticleVector)
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
            this->CollectStrayHits(pCluster, spanMinX, spanMaxX, collectedClusters);

	    if (!collectedClusters.empty())
	      this->AddInStrayClusters(pCluster, collectedClusters);
	    }
	}

    return (this->CreateThreeDParticles(protoParticleVector));
}

} // namespace lar_content

