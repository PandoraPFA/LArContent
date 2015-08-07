/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the end associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArObjects/LArPointingCluster.h"
#include "LArObjects/LArThreeDSlidingFitResult.h"

#include "LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef PfoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef PfoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

EndAssociatedPfosTool::EndAssociatedPfosTool() :
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexLongitudinalDistance(20.f),
    m_maxVertexTransverseDistance(3.5f),
    m_vertexAngularAllowance(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EndAssociatedPfosTool::Run(PfoHierarchyAlgorithm *const pAlgorithm, const Vertex *const /*pNeutrinoVertex*/, PfoInfoMap &pfoInfoMap)
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
            PfoInfoMap::iterator parentMapIter(pfoInfoMap.find(pParentPfo));

            if (pfoInfoMap.end() == parentMapIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            PfoInfo *const pParentPfoInfo(parentMapIter->second);
            const LArPointingCluster parentPointingCluster(*(pParentPfoInfo->GetSlidingFitResult3D()));
            const LArPointingCluster::Vertex &parentVertex(pParentPfoInfo->IsInnerLayerAssociated() ? parentPointingCluster.GetOuterVertex() : parentPointingCluster.GetInnerVertex());

            for (const ParticleFlowObject *const pPfo : unassignedPfos)
            {
                PfoInfoMap::iterator mapIter(pfoInfoMap.find(pPfo));

                if (pfoInfoMap.end() == mapIter)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                PfoInfo *const pPfoInfo(mapIter->second);

                const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));
                const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - parentVertex.GetPosition()).GetMagnitudeSquared() <
                    (pointingCluster.GetOuterVertex().GetPosition() - parentVertex.GetPosition()).GetMagnitudeSquared());

                const LArPointingCluster::Vertex &daughterVertex(useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());

                if (LArPointingClusterHelper::IsNode(parentVertex.GetPosition(), daughterVertex, m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                    LArPointingClusterHelper::IsNode(daughterVertex.GetPosition(), parentVertex, m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                    LArPointingClusterHelper::IsEmission(parentVertex.GetPosition(), daughterVertex,  m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
                    LArPointingClusterHelper::IsEmission(daughterVertex.GetPosition(), parentVertex,  m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance))
                {
                    associationsMade = true;
                    pParentPfoInfo->AddDaughterPfo(pPfoInfo->GetThisPfo());
                    pPfoInfo->SetParentPfo(pParentPfoInfo->GetThisPfo());
                    pPfoInfo->SetInnerLayerAssociation(useInner);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EndAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexAngularAllowance", m_vertexAngularAllowance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
