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
    m_nMaxTensorToolRepeats(1000),
    m_minClusterCaloHits(5),//5),
    m_searchRegion1D(5.f),    
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
}
        
//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    DeltaRayOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult));

    if (overlapResult.IsInitialized())
    {
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
        std::cout << "IN ALG SIZE: " << overlapResult.GetCommonMuonPfoList().size() << std::endl;
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
    
    std::cout << "IN FUNCTION LIST: " << commonMuonPfoList.size() << std::endl;
    std::cout << "IN FUNCTION MEMBER VARIABLE SIZE: " << overlapResult.GetCommonMuonPfoList().size() << std::endl;
    
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
    
    for (const ParticleFlowObject *const pNearbyMuonU : nearbyMuonPfosU)
    {
        for (const ParticleFlowObject *const pNearbyMuonV : nearbyMuonPfosV)
        {
            for (const ParticleFlowObject *const pNearbyMuonW : nearbyMuonPfosW)
            {
                if ((pNearbyMuonU == pNearbyMuonV) && (pNearbyMuonV == pNearbyMuonW))
                {
                    //std::cout << "FOUND COMMON MUON!" << std::endl;
                    //PandoraMonitoringApi::ViewEvent(this->GetPandora());                    
                    commonMuonPfoList.push_back(pNearbyMuonU);
                }
            }
        }
    }

    //std::cout << "HAVE NOT FOUND COMMON MUON!" << std::endl;
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());    
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
    
    return BaseAlgorithm::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
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

void ThreeViewDeltaRayMatchingAlgorithm::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pNewCluster));

    HitKDTree2D &kdTree((hitType == TPC_VIEW_U) ? m_kdTreeU : (hitType == TPC_VIEW_V) ? m_kdTreeV : m_kdTreeW);    
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);    

    CaloHitList caloHitList;
    pNewCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        hitToClusterMap[pCaloHit] = pNewCluster;
            
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
            ClusterList &invertedCloseClusters(clusterProximityMap.at(pNearbyCluster));
            
            if (std::find(invertedCloseClusters.begin(), invertedCloseClusters.end(), pNewCluster) == invertedCloseClusters.end())
                    invertedCloseClusters.push_back(pNewCluster);
        }
    }
    
    BaseAlgorithm::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pDeletedCluster));

    // Remove from hitToCluster map    
    HitToClusterMap &hitToClusterMap((hitType == TPC_VIEW_U) ? m_hitToClusterMapU : (hitType == TPC_VIEW_V) ? m_hitToClusterMapV : m_hitToClusterMapW);

    std::cout << "Begin deleting hits.. " << std::endl;
    
    const OrderedCaloHitList &orderedCaloHitList(pDeletedCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            auto iter(hitToClusterMap.find(pCaloHit));
            hitToClusterMap.erase(iter);
        }
    }

    std::cout << "Finished deleting hits.. " << std::endl;    

    // Remove from clusterProximity map
    ClusterProximityMap &clusterProximityMap((hitType == TPC_VIEW_U) ? m_clusterProximityMapU : (hitType == TPC_VIEW_V) ? m_clusterProximityMapV : m_clusterProximityMapW);

    const ClusterProximityMap::const_iterator clusterProximityIter(clusterProximityMap.find(pDeletedCluster));

    std::cout << "enter proximity... " << std::endl;    
    if (clusterProximityIter != clusterProximityMap.end())
    {
        ClusterVector closeClusterVector;

        for (const Cluster *const pCloseCluster : clusterProximityIter->second)
            closeClusterVector.push_back(pCloseCluster);

        std::sort(closeClusterVector.begin(), closeClusterVector.end(), LArClusterHelper::SortByNHits);

        std::cout << "closeClusterVector size: " << closeClusterVector.size() << std::endl;
        for (const Cluster *const pCloseCluster : closeClusterVector)
        {
            std::cout << "here" << std::endl;
            const ClusterProximityMap::iterator iter(clusterProximityMap.find(pCloseCluster));

            // May have also been deleted
            if (iter == clusterProximityMap.end())
                continue;
        
            ClusterList &invertedCloseClusters(iter->second);
            invertedCloseClusters.remove(pDeletedCluster);
        }

        clusterProximityMap.erase(clusterProximityIter);
    }

    std::cout << "finish proximity... " << std::endl;
    BaseAlgorithm::UpdateUponDeletion(pDeletedCluster);
}


} // namespace lar_content


/*

void ThreeViewDeltaRayMatchingAlgorithm::UpdateUponPfoCreation()
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
            pfoListName, pPfoList));

        if ((!pPfoList) || pPfoList->empty())
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            ClusterList pfoClusterList(pPfo->GetClusterList());

            for (const Cluster *const pCluster : pfoClusterList)
            {
                const HitType &hitType(LArClusterHelper::GetClusterHitType(pCluster));

                if ((hitType != TPC_VIEW_U) || (hitType != TPC_VIEW_V) || (hitType != TPC_VIEW_W))
                    continue;
                
                ClusterToPfoMap &clusterToPfoMap((hitType == TPC_VIEW_U) ? m_clusterToPfoMapU : (hitType == TPC_VIEW_V) ? m_clusterToPfoMapV : m_clusterToPfoMapW);

                if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
                    break;
                
                clusterToPfoMap[pCluster] = pPfo;
            }
        }
    }
}

*/
