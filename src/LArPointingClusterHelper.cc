/**
 *  @file   LArPointingClusterHelper.cc
 * 
 *  @brief  Implementation of the pointing cluster helper class.
 * 
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"
#include "Helpers/XmlHelper.h"

#include "LArPointingClusterHelper.h"

using namespace pandora;

namespace lar
{

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::IsNode(const CartesianVector &parentVertex, const CartesianVector &daughterVertex)
{
    if ((parentVertex - daughterVertex).GetMagnitudeSquared() < m_maxNodeRadiusSquared)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::IsEmission(const LArPointingCluster::Vertex &parentPointingVertex, const pandora::CartesianVector &daughterVertex)
{
    return LArPointingClusterHelper::IsPointing(daughterVertex, parentPointingVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::IsEmitted(const pandora::CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterPointingVertex)
{
    return LArPointingClusterHelper::IsPointing(parentVertex, daughterPointingVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPointingClusterHelper::IsPointing(const pandora::CartesianVector &vertex, const LArPointingCluster::Vertex &pointingVertex)
{
    const CartesianVector displacement(pointingVertex.GetPosition() - vertex);

    const float longitudinalDistance(pointingVertex.GetDirection().GetDotProduct(displacement));

    if ((longitudinalDistance < m_minPointingLongitudinalDistance) || (longitudinalDistance > m_maxPointingLongitudinalDistance))
        return false;

    const float maxTransverseDistanceSquared((m_maxPointingTransverseDistance * m_maxPointingTransverseDistance) + 
        (m_pointingAngularAllowance * m_pointingAngularAllowance * longitudinalDistance * longitudinalDistance));

    const float transverseDistanceSquared(pointingVertex.GetDirection().GetCrossProduct(displacement).GetMagnitudeSquared());

    if (transverseDistanceSquared > maxTransverseDistanceSquared)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPointingClusterHelper::m_maxNodeRadiusSquared = 2.5f * 2.5f;
float LArPointingClusterHelper::m_maxPointingLongitudinalDistance = 20.f;
float LArPointingClusterHelper::m_minPointingLongitudinalDistance = -2.5f;
float LArPointingClusterHelper::m_maxPointingTransverseDistance = 2.5f;
float LArPointingClusterHelper::m_pointingAngularAllowance = 0.0175f; // tan (1 degree)

StatusCode LArPointingClusterHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    float maxNodeRadius = 2.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxNodeRadius", maxNodeRadius));
    m_maxNodeRadiusSquared = maxNodeRadius * maxNodeRadius;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxPointingLongitudinalDistance", m_maxPointingLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinPointingLongitudinalDistance", m_minPointingLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxPointingTransverseDistance", m_maxPointingTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "PointingAngularAllowance", m_pointingAngularAllowance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
