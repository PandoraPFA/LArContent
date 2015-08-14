/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the branch associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArPointingCluster.h"
#include "LArObjects/LArThreeDSlidingFitResult.h"

#include "LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef NeutrinoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

BranchAssociatedPfosTool::BranchAssociatedPfosTool() :
    m_minNeutrinoVertexDistance(3.5),
    m_maxParentClusterDistance(3.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BranchAssociatedPfosTool::Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const Vertex *const pNeutrinoVertex, PfoInfoMap &pfoInfoMap)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    bool associationsMade(true);

    while (associationsMade)
    {
        associationsMade = false;
        PfoList assignedPfos, unassignedPfos;
        pAlgorithm->SeparatePfos(pfoInfoMap, assignedPfos, unassignedPfos);

        if (unassignedPfos.empty())
            break;

        for (const ParticleFlowObject *const pParentPfo : assignedPfos)
        {
            PfoInfo *const pParentPfoInfo(pfoInfoMap.at(pParentPfo));
            const Cluster *const pParentCluster3D(pParentPfoInfo->GetCluster3D());

            for (const ParticleFlowObject *const pPfo : unassignedPfos)
            {
                PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));
                const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));
                const float dNeutrinoVertex(LArClusterHelper::GetClosestDistance(pNeutrinoVertex->GetPosition(), pPfoInfo->GetCluster3D()));

                if (dNeutrinoVertex < m_minNeutrinoVertexDistance)
                    continue;

                const float dInnerVertex(LArClusterHelper::GetClosestDistance(pointingCluster.GetInnerVertex().GetPosition(), pParentCluster3D));
                const float dOuterVertex(LArClusterHelper::GetClosestDistance(pointingCluster.GetOuterVertex().GetPosition(), pParentCluster3D));

                if ((dInnerVertex < m_maxParentClusterDistance) || (dOuterVertex < m_maxParentClusterDistance))
                {
                    associationsMade = true;
                    pParentPfoInfo->AddDaughterPfo(pPfoInfo->GetThisPfo());
                    pPfoInfo->SetParentPfo(pParentPfoInfo->GetThisPfo());
                    pPfoInfo->SetInnerLayerAssociation(dInnerVertex < dOuterVertex);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNeutrinoVertexDistance", m_minNeutrinoVertexDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxParentClusterDistance", m_maxParentClusterDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
