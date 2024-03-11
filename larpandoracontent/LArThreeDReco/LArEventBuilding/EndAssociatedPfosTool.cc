/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.cc
 *
 *  @brief  Implementation of the end associated pfos tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef NeutrinoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

EndAssociatedPfosTool::EndAssociatedPfosTool() :
    m_minNeutrinoVertexDistance(5.f),
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexLongitudinalDistance(20.f),
    m_maxVertexTransverseDistance(3.5f),
    m_vertexAngularAllowance(3.f),
    m_maxParentEndpointDistance(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EndAssociatedPfosTool::Run(const NeutrinoHierarchyAlgorithm *const pAlgorithm, const Vertex *const pNeutrinoVertex, PfoInfoMap &pfoInfoMap)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool associationsMade(true);

    while (associationsMade)
    {
        associationsMade = false;
        PfoVector assignedPfos, unassignedPfos;
        pAlgorithm->SeparatePfos(pfoInfoMap, assignedPfos, unassignedPfos);

        if (unassignedPfos.empty())
            break;

        // ATTN May want to reconsider precise association mechanics for complex situations
        PfoSet recentlyAssigned;

        for (const ParticleFlowObject *const pParentPfo : assignedPfos)
        {
            PfoInfo *const pParentPfoInfo(pfoInfoMap.at(pParentPfo));
            const LArPointingCluster parentPointingCluster(*(pParentPfoInfo->GetSlidingFitResult3D()));

            const LArPointingCluster::Vertex &parentEndpoint(
                pParentPfoInfo->IsInnerLayerAssociated() ? parentPointingCluster.GetOuterVertex() : parentPointingCluster.GetInnerVertex());
            const float neutrinoVertexDistance((parentEndpoint.GetPosition() - pNeutrinoVertex->GetPosition()).GetMagnitude());

            if (neutrinoVertexDistance < m_minNeutrinoVertexDistance)
                continue;

            for (const ParticleFlowObject *const pPfo : unassignedPfos)
            {
                if (recentlyAssigned.count(pPfo))
                    continue;

                PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

                const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));
                const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - parentEndpoint.GetPosition()).GetMagnitudeSquared() <
                    (pointingCluster.GetOuterVertex().GetPosition() - parentEndpoint.GetPosition()).GetMagnitudeSquared());

                const LArPointingCluster::Vertex &daughterVertex(useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());

                if (LArPointingClusterHelper::IsNode(parentEndpoint.GetPosition(), daughterVertex, m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                    LArPointingClusterHelper::IsNode(daughterVertex.GetPosition(), parentEndpoint, m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                    LArPointingClusterHelper::IsEmission(parentEndpoint.GetPosition(), daughterVertex, m_minVertexLongitudinalDistance,
                        m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
                    LArPointingClusterHelper::IsEmission(daughterVertex.GetPosition(), parentEndpoint, m_minVertexLongitudinalDistance,
                        m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
                    this->IsCloseToParentEndpoint(parentEndpoint.GetPosition(), pParentPfoInfo->GetCluster3D(), pPfoInfo->GetCluster3D()))
                {
                    associationsMade = true;
                    pParentPfoInfo->AddDaughterPfo(pPfoInfo->GetThisPfo());
                    pPfoInfo->SetParentPfo(pParentPfoInfo->GetThisPfo());
                    pPfoInfo->SetInnerLayerAssociation(useInner);
                    recentlyAssigned.insert(pPfoInfo->GetThisPfo());
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EndAssociatedPfosTool::IsCloseToParentEndpoint(
    const CartesianVector &parentEndpoint, const Cluster *const pParentCluster3D, const Cluster *const pDaughterCluster3D) const
{
    try
    {
        CartesianVector parentPosition3D(0.f, 0.f, 0.f), daughterPosition3D(0.f, 0.f, 0.f);
        LArClusterHelper::GetClosestPositions(pParentCluster3D, pDaughterCluster3D, parentPosition3D, daughterPosition3D);

        if (((parentPosition3D - parentEndpoint).GetMagnitude() < m_maxParentEndpointDistance) &&
            ((parentPosition3D - daughterPosition3D).GetMagnitude() < m_maxParentEndpointDistance))
            return true;
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EndAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinNeutrinoVertexDistance", m_minNeutrinoVertexDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexAngularAllowance", m_vertexAngularAllowance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxParentEndpointDistance", m_maxParentEndpointDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
