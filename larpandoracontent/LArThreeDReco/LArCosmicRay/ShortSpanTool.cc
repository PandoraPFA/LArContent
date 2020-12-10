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
    m_minXOverlapFraction(0.7f)
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

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::InvestigateShortSpans(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
    for (TensorType::Element &element : elementList)
    {
        const Cluster *pShortCluster(nullptr);

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
	PandoraMonitoringApi::ViewEvent(this->GetPandora());
	*/
        ////////////////////////// 

        // Check if short span creates ambiguities
        if(!this->GetShortCluster(element, pShortCluster))
            continue;

        HitType badHitType(LArClusterHelper::GetClusterHitType(pShortCluster));

        // Create list of stray clusters - those that do not pass tensor thresholds
        pAlgorithm->InitialiseStrayClusterList(badHitType);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &strayClusterList, "STRAY CLUSTERS", BLACK);
        
        // Collect any clusters to add into the short cluster 
	float spanMinX(0.f), spanMaxX(0.f);
	this->GetGoodXOverlapExtrema(element, badHitType, spanMinX, spanMaxX);

	ClusterList collectedClusters;
	pAlgorithm->CollectStrayHits(pShortCluster, spanMinX, spanMaxX, collectedClusters);
        
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &collectedClusters, "COLLECTED CLUSTERS", DARKGREEN);

        // Add stray clusters in
        if (!collectedClusters.empty())
        {
	  pAlgorithm->AddInStrayClusters(pShortCluster, collectedClusters);
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








