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
    m_maxSeparation(2.f),
    m_minNConnectedClusters(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UnambiguousDeltaRayTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    m_pParentAlgorithm = pAlgorithm;

    if (PandoraContentApi::GetSettings(*m_pParentAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);

    return this->ExamineUnambiguousElements(elementList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UnambiguousDeltaRayTool::ExamineUnambiguousElements(TensorType::ElementList &elementList)
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::Element &element : elementList)
    {
        if (!this->IsConnected(element))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.emplace_back(element.GetClusterU());
        protoParticle.m_clusterList.emplace_back(element.GetClusterV());
        protoParticle.m_clusterList.emplace_back(element.GetClusterW());
        protoParticleVector.emplace_back(protoParticle);
    }

    return m_pParentAlgorithm->CreatePfos(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UnambiguousDeltaRayTool::IsConnected(const TensorType::Element &element) const
{
    PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
    {
        unsigned int connectedClusterCount(0);

        for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            ClusterList muonClusterList;
            LArPfoHelper::GetClusters(pMuonPfo, hitType, muonClusterList);

            const float separation(LArClusterHelper::GetClosestDistance(element.GetCluster(hitType), muonClusterList));

            if (separation < m_maxSeparation)
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinNConnectedClusters", m_minNConnectedClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
