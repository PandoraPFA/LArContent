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

bool SpanAcceptanceTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &/*overlapTensor*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    TensorType::ElementList elementList;
    pAlgorithm->GetUnambiguousElements(true, elementList);

    if (elementList.empty())
        return changesMade;

    this->InvestigateSpanAcceptances(pAlgorithm, elementList, changesMade);

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SpanAcceptanceTool::InvestigateSpanAcceptances(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
{
  ProtoParticleVector protoParticleVector;
    for (TensorType::Element &element : elementList)
    {
        if (!this->IsConnected(element))
            continue;

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

    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
    {
        unsigned int connectedClusterCount(0);
        for (const HitType &hitType : hitTypeVector)
        {
            ClusterList muonClusterList;
            LArPfoHelper::GetClusters(pMuonPfo, hitType, muonClusterList);

            if (muonClusterList.size() != 1)
            {
                std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
                continue;
            }

            const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList));

            if (separation < 2.f)
                ++connectedClusterCount;
        }

        if (connectedClusterCount > 1)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode SpanAcceptanceTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content








