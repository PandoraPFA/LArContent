/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/UnambiguousDeltaRayTool.cc
 *
 *  @brief  Implementation of the unambiguous delta ray tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/UnambiguousDeltaRayTool.h"

using namespace pandora;

namespace lar_content
{

UnambiguousDeltaRayTool::UnambiguousDeltaRayTool() :
    m_minSeparation(2.f),
    m_minNConnectedClusters(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UnambiguousDeltaRayTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);

    this->ExamineUnambiguousElements(pAlgorithm, elementList, changesMade);

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UnambiguousDeltaRayTool::ExamineUnambiguousElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType::ElementList &elementList, bool &changesMade)
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

bool UnambiguousDeltaRayTool::IsConnected(const TensorType::Element &element) const
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
    {
        unsigned int connectedClusterCount(0);
        
        for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            ClusterList muonClusterList;
            LArPfoHelper::GetClusters(pMuonPfo, hitType, muonClusterList);

            const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList));

            if (separation < m_minSeparation)
                ++connectedClusterCount;
        }

        if (connectedClusterCount > m_minNConnectedClusters)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode UnambiguousDeltaRayTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparation", m_minSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNConnectedClusters", m_minNConnectedClusters));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content








