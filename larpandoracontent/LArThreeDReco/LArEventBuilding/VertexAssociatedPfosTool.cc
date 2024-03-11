/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.cc
 *
 *  @brief  Implementation of the vertex associated pfos tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef NeutrinoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

VertexAssociatedPfosTool::VertexAssociatedPfosTool() :
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexLongitudinalDistance(20.f),
    m_maxVertexTransverseDistance(3.5f),
    m_vertexAngularAllowance(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexAssociatedPfosTool::Run(const NeutrinoHierarchyAlgorithm *const pAlgorithm, const Vertex *const pNeutrinoVertex, PfoInfoMap &pfoInfoMap)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const CartesianVector &neutrinoVertex(pNeutrinoVertex->GetPosition());

    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoInfoMap)
        sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo : sortedPfos)
    {
        PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

        if (pPfoInfo->IsNeutrinoVertexAssociated() || pPfoInfo->GetParentPfo())
            continue;

        const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));
        const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - neutrinoVertex).GetMagnitudeSquared() <
            (pointingCluster.GetOuterVertex().GetPosition() - neutrinoVertex).GetMagnitudeSquared());

        const LArPointingCluster::Vertex &daughterVertex(useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());

        if (LArPointingClusterHelper::IsNode(neutrinoVertex, daughterVertex, m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
            LArPointingClusterHelper::IsEmission(neutrinoVertex, daughterVertex, m_minVertexLongitudinalDistance,
                m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance))
        {
            pPfoInfo->SetNeutrinoVertexAssociation(true);
            pPfoInfo->SetInnerLayerAssociation(useInner);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexAngularAllowance", m_vertexAngularAllowance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
