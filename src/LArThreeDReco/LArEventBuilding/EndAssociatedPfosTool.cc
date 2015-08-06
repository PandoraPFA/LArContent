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
    m_maxVertexTransverseDistance(1.5f),
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
        this->SeparatePfos(pfoInfoMap, assignedPfos, unassignedPfos);

        if (unassignedPfos.empty())
            break;

        for (const ParticleFlowObject *const pParentPfo : assignedPfos)
        {
            PfoInfoMap::iterator parentMapIter(pfoInfoMap.find(pParentPfo));

            if (pfoInfoMap.end() == parentMapIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            PfoInfo *const pParentPfoInfo(parentMapIter->second);
            const LArPointingCluster parentPointingCluster(*(pParentPfoInfo->GetSlidingFitResult3D()));
            const LArPointingCluster::Vertex &parentVertex(pParentPfoInfo->IsInnerLayerAssociated() ? parentPointingCluster.GetInnerVertex() : parentPointingCluster.GetOuterVertex());

            for (const ParticleFlowObject *const pDaughterPfo : unassignedPfos)
            {
                PfoInfoMap::iterator daughterMapIter(pfoInfoMap.find(pDaughterPfo));

                if (pfoInfoMap.end() == daughterMapIter)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                PfoInfo *const pDaughterPfoInfo(daughterMapIter->second);
                // TODO, proceed as for vertex associated pfos
                std::cout << " pDaughterPfoInfo " << pDaughterPfoInfo << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EndAssociatedPfosTool::SeparatePfos(const PfoInfoMap &pfoInfoMap, PfoList &assignedPfos, PfoList &unassignedPfos) const
{
    for (const PfoInfoMap::value_type mapIter : pfoInfoMap)
    {
        const PfoInfo *const pPfoInfo(mapIter.second);

        if (pPfoInfo->IsNeutrinoVertexAssociated() || pPfoInfo->GetParentPfo())
        {
            assignedPfos.insert(pPfoInfo->GetThisPfo());
        }
        else
        {
            unassignedPfos.insert(pPfoInfo->GetThisPfo());
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
