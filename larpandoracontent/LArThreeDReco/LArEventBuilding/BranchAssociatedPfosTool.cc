/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.cc
 *
 *  @brief  Implementation of the branch associated pfos tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef NeutrinoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

BranchAssociatedPfosTool::BranchAssociatedPfosTool() :
    m_minNeutrinoVertexDistance(5.f), m_trackBranchAdditionFraction(0.4f), m_maxParentClusterDistance(3.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BranchAssociatedPfosTool::Run(const NeutrinoHierarchyAlgorithm *const pAlgorithm, const Vertex *const pNeutrinoVertex, PfoInfoMap &pfoInfoMap)
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
            const Cluster *const pParentCluster3D(pParentPfoInfo->GetCluster3D());

            const bool parentIsTrack(LArPfoHelper::IsTrack(pParentPfo));
            const ThreeDSlidingFitResult &parentFitResult(*(pParentPfoInfo->GetSlidingFitResult3D()));

            const float parentLength3D((parentFitResult.GetGlobalMinLayerPosition() - parentFitResult.GetGlobalMaxLayerPosition()).GetMagnitude());
            const CartesianVector &parentVertexPosition(pParentPfoInfo->IsInnerLayerAssociated() ? parentFitResult.GetGlobalMinLayerPosition()
                                                                                                 : parentFitResult.GetGlobalMaxLayerPosition());

            for (const ParticleFlowObject *const pPfo : unassignedPfos)
            {
                if (recentlyAssigned.count(pPfo))
                    continue;

                PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));
                const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));

                const float dNeutrinoVertex(std::min((pointingCluster.GetInnerVertex().GetPosition() - pNeutrinoVertex->GetPosition()).GetMagnitude(),
                    (pointingCluster.GetOuterVertex().GetPosition() - pNeutrinoVertex->GetPosition()).GetMagnitude()));

                if (dNeutrinoVertex < m_minNeutrinoVertexDistance)
                    continue;

                const float dParentVertex(std::min((pointingCluster.GetInnerVertex().GetPosition() - parentVertexPosition).GetMagnitude(),
                    (pointingCluster.GetOuterVertex().GetPosition() - parentVertexPosition).GetMagnitude()));

                if (parentIsTrack && (dParentVertex < m_trackBranchAdditionFraction * parentLength3D))
                    continue;

                const float dInnerVertex(LArClusterHelper::GetClosestDistance(pointingCluster.GetInnerVertex().GetPosition(), pParentCluster3D));
                const float dOuterVertex(LArClusterHelper::GetClosestDistance(pointingCluster.GetOuterVertex().GetPosition(), pParentCluster3D));

                if ((dInnerVertex < m_maxParentClusterDistance) || (dOuterVertex < m_maxParentClusterDistance))
                {
                    associationsMade = true;
                    pParentPfoInfo->AddDaughterPfo(pPfoInfo->GetThisPfo());
                    pPfoInfo->SetParentPfo(pParentPfoInfo->GetThisPfo());
                    pPfoInfo->SetInnerLayerAssociation(dInnerVertex < dOuterVertex);
                    recentlyAssigned.insert(pPfoInfo->GetThisPfo());
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinNeutrinoVertexDistance", m_minNeutrinoVertexDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrackBranchAdditionFraction", m_trackBranchAdditionFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxParentClusterDistance", m_maxParentClusterDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
