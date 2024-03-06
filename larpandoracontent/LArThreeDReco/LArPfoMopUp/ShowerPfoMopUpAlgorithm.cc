/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerPfoMopUpAlgorithm.cc
 *
 *  @brief  Implementation of the shower pfo mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerPfoMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerPfoMopUpAlgorithm::ShowerPfoMopUpAlgorithm() :
    VertexBasedPfoMopUpAlgorithm(), m_maxVertexLongitudinalDistance(20.f), m_vertexAngularAllowance(3.f)
{
    // ATTN Some default values differ from base class
    m_maxVertexTransverseDistance = 3.5f;
    m_meanBoundedFractionCut = 0.5f;
    m_maxBoundedFractionCut = 0.7f;
    m_minBoundedFractionCut = 0.f;
    m_minConsistentDirections = 1;
    m_minConsistentDirectionsTrack = 2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerPfoMopUpAlgorithm::IsVertexAssociated(const CartesianVector &vertex2D, const LArPointingCluster &pointingCluster) const
{
    return (LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
            LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
            LArPointingClusterHelper::IsEmission(vertex2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance,
                m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
            LArPointingClusterHelper::IsEmission(vertex2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance,
                m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerPfoMopUpAlgorithm::PfoAssociation ShowerPfoMopUpAlgorithm::GetPfoAssociation(
    const Pfo *const pVertexPfo, const Pfo *const pDaughterPfo, HitTypeToAssociationMap &hitTypeToAssociationMap) const
{
    return PfoAssociation(pVertexPfo, pDaughterPfo, hitTypeToAssociationMap[TPC_VIEW_U], hitTypeToAssociationMap[TPC_VIEW_V],
        hitTypeToAssociationMap[TPC_VIEW_W]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerPfoMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexAngularAllowance", m_vertexAngularAllowance));

    return VertexBasedPfoMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
