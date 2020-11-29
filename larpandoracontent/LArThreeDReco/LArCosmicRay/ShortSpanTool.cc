/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ShortSpanTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ShortSpanTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

ShortSpanTool::ShortSpanTool() :
    m_minXOverlapFraction(0.7f),
    m_isStrayListUInitialised(false),
    m_isStrayListVInitialised(false),
    m_isStrayListWInitialised(false),
    m_strayClusterListU(ClusterList()),
    m_strayClusterListV(ClusterList()),
    m_strayClusterListW(ClusterList())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

    bool ShortSpanTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->InvestigateShortSpans(pAlgorithm, elementList, changesMade);

    this->ClearStrayClusterLists();
    
    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::InvestigateShortSpans(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
    for (TensorType::Element &element : elementList)
    {
        const Cluster *pShortCluster(nullptr);

        // Check if short span creates ambiguities
        if(!this->GetShortCluster(element, pShortCluster))
            return;

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
        */
        ////////////////////////// 

        // Create list of stray clusters - those that do not pass tensor thresholds
        this->InitialiseStrayClusterList(pAlgorithm, LArClusterHelper::GetClusterHitType(pShortCluster));

        HitType hitType(LArClusterHelper::GetClusterHitType(pShortCluster));
        ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &strayClusterList, "STRAY CLUSTERS", BLACK);
        
        // Collect any clusters to add into the short cluster 
        ClusterList collectedClusters;
        this->CollectStrayHits(element, pShortCluster, strayClusterList, collectedClusters);
        
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &collectedClusters, "COLLECTED CLUSTERS", DARKGREEN);

        // Add stray clusters in
        if (!collectedClusters.empty())
        {
            pAlgorithm->UpdateUponDeletion(pShortCluster);

            for(const Cluster *const pCollectedCluster : collectedClusters)
            {
                const ClusterList::const_iterator deleteIter(std::find(strayClusterList.begin(), strayClusterList.end(), pCollectedCluster));
                strayClusterList.erase(deleteIter);
                
                pAlgorithm->UpdateUponDeletion(pCollectedCluster);

                std::string clusterListName(pAlgorithm->GetClusterListName(LArClusterHelper::GetClusterHitType(pShortCluster)));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pShortCluster, pCollectedCluster, clusterListName, clusterListName));
            }

            pAlgorithm->UpdateForNewCluster(pShortCluster);

            changesMade = true;
        }

        //////////////////////////
        /*        
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &uCluster, "uCluster_2", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &vCluster, "vCluster_2", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &wCluster, "wCluster_2", VIOLET);        

        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        */
        //////////////////////////
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::InitialiseStrayClusterList(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const HitType &hitType)
{
    if (this->IsStrayClusterListInitialised(hitType))
        return;
    
    ClusterList &strayClusterList((hitType == TPC_VIEW_U) ? m_strayClusterListU : (hitType == TPC_VIEW_V) ? m_strayClusterListV : m_strayClusterListW);
    
    for (const Cluster *const pCluster : pAlgorithm->GetInputClusterList(hitType))
    {
        if (!pAlgorithm->DoesClusterPassTesorThreshold(pCluster))
            strayClusterList.push_back(pCluster);
    }

    bool &isStrayClusterListInitialised((hitType == TPC_VIEW_U) ? m_isStrayListUInitialised : (hitType == TPC_VIEW_V) ? m_isStrayListVInitialised : m_isStrayListWInitialised);
    isStrayClusterListInitialised = true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

bool ShortSpanTool::GetShortCluster(const TensorType::Element &element, const Cluster *&pShortCluster) const
{
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    float shortestSpan(std::numeric_limits<float>::max()), longestSpan(-std::numeric_limits<float>::max());
    for (const HitType &hitType : hitTypeVector)
    {
        const float viewXSpan(element.GetOverlapResult().GetViewXSpan(hitType));
        
        if (viewXSpan < shortestSpan)
        {
            shortestSpan = viewXSpan;
            pShortCluster = element.GetCluster(hitType);
        }

        if(viewXSpan > longestSpan)
            longestSpan = viewXSpan;
    }

    // Check if failed because of span inconsistencies
    if ((longestSpan < std::numeric_limits<float>::epsilon()) || (element.GetOverlapResult().GetXOverlap().GetXOverlapSpan() / longestSpan > m_minXOverlapFraction))
        return false;

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

    // Was the failure a result of a short span?
    return ((middleSpan - shortestSpan) > (longestSpan - middleSpan));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::CollectStrayHits(const TensorType::Element &element, const Cluster *const pShortCluster, const ClusterList &strayClusterList, ClusterList &collectedClusters) const
{
    HitType badHitType(LArClusterHelper::GetClusterHitType(pShortCluster));
    
    float spanMinX(0.f), spanMaxX(0.f);
    this->GetGoodXOverlapExtrema(element, badHitType, spanMinX, spanMaxX);
    
    for (const Cluster *const pCluster : strayClusterList)
    {
        float xMin(-std::numeric_limits<float>::max()), xMax(+std::numeric_limits<float>::max());
        pCluster->GetClusterSpanX(xMin, xMax);

        if (!(((xMin > spanMinX) && (xMin < spanMaxX)) || ((xMax > spanMinX) && (xMax < spanMaxX))))
            continue;
        
        if (LArClusterHelper::GetClosestDistance(pShortCluster, pCluster) > 2.f)
            continue;

        collectedClusters.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::ClearStrayClusterLists()
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

bool ShortSpanTool::IsStrayClusterListInitialised(const HitType &hitType) const
{
    if ((hitType != TPC_VIEW_U) && (hitType != TPC_VIEW_V) && (hitType != TPC_VIEW_W))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    return (hitType == TPC_VIEW_U) ? m_isStrayListUInitialised : (hitType == TPC_VIEW_V) ? m_isStrayListVInitialised : m_isStrayListWInitialised;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::GetGoodXOverlapExtrema(const TensorType::Element &element, const HitType &badHitType, float &minX, float &maxX) const
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
    
StatusCode ShortSpanTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
















//------------------------------------------------------------------------------------------------------------------------------------------
/*
void ShortSpanTool::CollectMuonHits(const TensorType::Element &element, const Cluster *const pShortCluster, const ClusterList &strayClusterList, ClusterList &collectedClusters) const
{
    HitType badHitType(LArClusterHelper::GetClusterHitType(pShortCluster));
    
    float spanMinX(0.f), spanMaxX(0.f);
    this->GetGoodXOverlapExtrema(element, badHitType, spanMinX, spanMaxX);

    const PfoList &commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
    {
        std::cout << "ISOBEL: MORE THAN ONE COMMON MUON IN COMMON MUON LIST" << std::endl;
        return;
    }

    ClusterList muonCluster;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), badHitType, muonCluster);

    for (const Cluster *const pCluster : strayClusterList)
    {
        float xMin(-std::numeric_limits<float>::max()), xMax(+std::numeric_limits<float>::max());
        pCluster->GetClusterSpanX(xMin, xMax);

        if (!(((xMin > spanMinX) && (xMin < spanMaxX)) || ((xMax > spanMinX) && (xMax < spanMaxX))))
            continue;
        
        if (LArClusterHelper::GetClosestDistance(pShortCluster, pCluster) > 2.f)
            continue;

        collectedClusters.push_back(pCluster);
    }
}
*/
