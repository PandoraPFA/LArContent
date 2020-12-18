/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRayg/SpandAcceptanceTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/SpanAcceptanceTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

SpanAcceptanceTool::SpanAcceptanceTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SpanAcceptanceTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->InvestigateSpanAcceptances(pAlgorithm, elementList, changesMade);

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SpanAcceptanceTool::InvestigateSpanAcceptances(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
  ProtoParticleVector protoParticleVector;
    for (TensorType::Element &element : elementList)
    {
        //////////////////////////
      
        std::cout << "uSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_U) << std::endl;
        std::cout << "vSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_V) << std::endl;
        std::cout << "wSpan: " << element.GetOverlapResult().GetViewXSpan(TPC_VIEW_W) << std::endl;
        std::cout << "overlap: " << element.GetOverlapResult().GetXOverlap().GetXOverlapSpan() << std::endl;
                
        ClusterList uCluster({element.GetCluster(TPC_VIEW_U)}), vCluster({element.GetCluster(TPC_VIEW_V)}), wCluster({element.GetCluster(TPC_VIEW_W)});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &uCluster, "uCluster_1", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &vCluster, "vCluster_1", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &wCluster, "wCluster_1", VIOLET);

        ////////////////////////// 

	// Is there one common muon?
	PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

	if (commonMuonPfoList.size() != 1)
        {
	  std::cout << "too many common muons" << std::endl;
	  PandoraMonitoringApi::ViewEvent(this->GetPandora());
	  continue;
	}

	if (!this->IsConnected(element))
        {
	  std::cout << "not enougth connection points" << std::endl;
	  PandoraMonitoringApi::ViewEvent(this->GetPandora());
	  continue;
	}
	PandoraMonitoringApi::ViewEvent(this->GetPandora());


        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(element.GetClusterU());
        protoParticle.m_clusterList.push_back(element.GetClusterV());
        protoParticle.m_clusterList.push_back(element.GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    changesMade |= pAlgorithm->CreatePfos(protoParticleVector);

}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool SpanAcceptanceTool::IsConnected(const TensorType::Element &element) const
{
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

     unsigned int connectedClusterCount(0);
    for (const HitType &hitType : hitTypeVector)
    {
        ClusterList muonClusterList;
        LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);

        if (muonClusterList.size() != 1)
        {
            std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
            continue;
        }

        const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList));

        if (separation < 2.f)
            ++connectedClusterCount;
    }

    return (connectedClusterCount > 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode SpanAcceptanceTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content








