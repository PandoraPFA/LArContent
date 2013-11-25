/**
 *  @file   LArContent/src/LArHelpers/LArVertexHelper.cc
 * 
 *  @brief  Implementation of the vertex helper class.
 * 
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"
#include "Helpers/XmlHelper.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArThreeDHelper.h"
#include "LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar
{

LArVertexHelper::NameToVertexMap LArVertexHelper::m_nameToVertexMap;
std::string LArVertexHelper::m_currentVertexName;

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::AddVertex(const std::string &vertexName, const CartesianVector &vertex)
{
    if (!(m_nameToVertexMap.insert(NameToVertexMap::value_type(vertexName, vertex)).second))
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::RemoveVertex(const std::string &vertexName)
{
    NameToVertexMap::iterator iter = m_nameToVertexMap.find(vertexName);

    if (m_nameToVertexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    m_nameToVertexMap.erase(iter);

    if (vertexName == m_currentVertexName)
        m_currentVertexName.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::UpdateVertex(const std::string &vertexName, const CartesianVector &vertex)
{
    NameToVertexMap::iterator iter = m_nameToVertexMap.find(vertexName);

    if (m_nameToVertexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    iter->second = vertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::SetCurrentVertex(const std::string &vertexName)
{
    NameToVertexMap::const_iterator iter = m_nameToVertexMap.find(vertexName);

    if (m_nameToVertexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    m_currentVertexName = vertexName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArVertexHelper::GetVertex(const std::string &vertexName)
{
    NameToVertexMap::const_iterator iter = m_nameToVertexMap.find(vertexName);

    if (m_nameToVertexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArVertexHelper::GetCurrentVertexName()
{
    NameToVertexMap::const_iterator iter = m_nameToVertexMap.find(m_currentVertexName);

    if (m_nameToVertexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_currentVertexName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &LArVertexHelper::GetCurrentVertex()
{
    NameToVertexMap::const_iterator iter = m_nameToVertexMap.find(m_currentVertexName);

    if (m_nameToVertexMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsForwardInZ(const Cluster *const pCluster, const CartesianVector &specifiedVertex)
{
    CartesianVector startPosition(0.f, 0.f, 0.f), startDirection(0.f, 0.f, 0.f);
    LArVertexHelper::GetInnerVertexAndDirection(pCluster, startPosition, startDirection); 

    const CartesianVector endPosition(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));
    return LArVertexHelper::IsDirectionCorrect(specifiedVertex, startPosition, endPosition, startDirection);    
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsBackwardInZ(const Cluster *const pCluster, const CartesianVector &specifiedVertex)
{
    CartesianVector startPosition(0.f, 0.f, 0.f), startDirection(0.f, 0.f, 0.f);
    LArVertexHelper::GetOuterVertexAndDirection(pCluster, startPosition, startDirection);

    const CartesianVector endPosition(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
    return LArVertexHelper::IsDirectionCorrect(specifiedVertex, startPosition, endPosition, startDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsBackwardInZ3D(const Cluster *const pCluster)
{
    const CartesianVector theVertex2D(LArGeometryHelper::ProjectPosition(LArVertexHelper::GetCurrentVertex(), LArThreeDHelper::GetClusterHitType(pCluster)));
    return LArVertexHelper::IsBackwardInZ(pCluster, theVertex2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsForwardInZ3D(const Cluster *const pCluster)
{
    const CartesianVector theVertex2D(LArGeometryHelper::ProjectPosition(LArVertexHelper::GetCurrentVertex(), LArThreeDHelper::GetClusterHitType(pCluster)));
    return LArVertexHelper::IsForwardInZ(pCluster, theVertex2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterDirection LArVertexHelper::GetDirectionInZ(const Cluster *const pCluster)
{
    if (!LArVertexHelper::DoesCurrentVertexExist())
        return DIRECTION_UNKNOWN;

    // calculate inner and outer position and direction
    CartesianVector innerPosition(0.f, 0.f, 0.f), innerDirection(0.f, 0.f, 0.f);
    LArVertexHelper::GetInnerVertexAndDirection(pCluster, innerPosition, innerDirection);

    CartesianVector outerPosition(0.f, 0.f, 0.f), outerDirection(0.f, 0.f, 0.f);
    LArVertexHelper::GetOuterVertexAndDirection(pCluster, outerPosition, outerDirection);

    // calculate direction using closest end to vertex
    bool isForwardInZ(false), isBackwardInZ(false);
    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex()); 

    if ((innerPosition - theVertex).GetMagnitudeSquared() < (outerPosition - theVertex).GetMagnitudeSquared())
    {
        isForwardInZ = LArVertexHelper::IsDirectionCorrect(LArVertexHelper::GetCurrentVertex(), innerPosition, outerPosition, innerDirection);
    }
    else
    {
        isBackwardInZ = LArVertexHelper::IsDirectionCorrect(LArVertexHelper::GetCurrentVertex(), outerPosition, innerPosition, outerDirection);
    }

    // analyse the results
    if (isForwardInZ && !isBackwardInZ)
        return FORWARD;

    if (!isForwardInZ && isBackwardInZ)
        return BACKWARD;

    return DIRECTION_AMBIGUOUS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArVertexHelper::GetDistanceSquaredToCurrentVertex(const Cluster *const pCluster)
{
    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());

    CartesianVector innerPosition(0.f, 0.f, 0.f), innerDirection(0.f, 0.f, 0.f);
    CartesianVector outerPosition(0.f, 0.f, 0.f), outerDirection(0.f, 0.f, 0.f);

    LArVertexHelper::GetInnerVertexAndDirection(pCluster, innerPosition, innerDirection); 
    LArVertexHelper::GetOuterVertexAndDirection(pCluster, outerPosition, outerDirection); 

    const float innerDistanceSquared((innerPosition-theVertex).GetMagnitudeSquared());
    const float outerDistanceSquared((outerPosition-theVertex).GetMagnitudeSquared());

    return std::min(innerDistanceSquared, outerDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArVertexHelper::GetImpactParameterToCurrentVertex(const Cluster *const pCluster)
{
    // Calculate impact parameters
    float longitudinal(-std::numeric_limits<float>::max());
    float transverse(+std::numeric_limits<float>::max());

    LArVertexHelper::GetImpactParametersToCurrentVertex(pCluster, longitudinal, transverse);

    // Vertex must be upstream of cluster
    if ( longitudinal > -m_longitudinalWindow ) 
        return transverse;
        
    // Cluster doesn't point back to vertex
    return std::numeric_limits<float>::max();;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::GetImpactParametersToCurrentVertex(const Cluster *const pCluster, float &longitudinal, float &transverse)
{
    // Calculate inner and outer vertex and direction
    CartesianVector innerPosition(0.f, 0.f, 0.f), innerDirection(0.f, 0.f, 0.f);
    CartesianVector outerPosition(0.f, 0.f, 0.f), outerDirection(0.f, 0.f, 0.f);

    LArVertexHelper::GetInnerVertexAndDirection(pCluster, innerPosition, innerDirection); 
    LArVertexHelper::GetOuterVertexAndDirection(pCluster, outerPosition, outerDirection); 

    // Calculate impact parameters using closest end to vertex
    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex()); 

    if ((innerPosition - theVertex).GetMagnitudeSquared() < (outerPosition - theVertex).GetMagnitudeSquared())
    {
        return LArVertexHelper::GetImpactParameters(innerPosition, innerDirection, theVertex, longitudinal, transverse);
    }
    else
    {
        return LArVertexHelper::GetImpactParameters(outerPosition, outerDirection, theVertex, longitudinal, transverse);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsDirectionCorrect(const CartesianVector &specifiedVertex, const CartesianVector &startPosition, const CartesianVector &endPosition,
    const CartesianVector &direction)
{
    float rL(0.f), rT(0.f);
    LArVertexHelper::GetImpactParameters(startPosition, direction, specifiedVertex, rL, rT);
    const float r((startPosition - endPosition).GetMagnitude());

    // Vertex must be inside a 60 degree triangle with apex placed one third of the way along the cluster (i.e. rL > +rT / tan(60) - r / 3.)
    if (rL > +0.57735f * rT - 0.33333f * r)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::GetInnerVertexAndDirection(const Cluster *const pCluster, CartesianVector &innerPosition, CartesianVector &innerDirection)
{
    ClusterHelper::ClusterFitResult innerLayerFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pCluster, m_numberOfLayersToFit, innerLayerFit));

    const CartesianVector &innerIntercept(innerLayerFit.GetIntercept());
    const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));

    innerDirection = innerLayerFit.GetDirection();
    innerPosition = innerIntercept + innerDirection * (innerDirection.GetDotProduct(innerCentroid - innerIntercept));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::GetOuterVertexAndDirection(const Cluster *const pCluster, CartesianVector &outerPosition, CartesianVector &outerDirection)
{
    ClusterHelper::ClusterFitResult outerLayerFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pCluster, m_numberOfLayersToFit, outerLayerFit));

    const CartesianVector &outerIntercept(outerLayerFit.GetIntercept());
    const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

    outerDirection = outerLayerFit.GetDirection() * -1.f; // Note: minus sign!
    outerPosition  = outerIntercept + outerDirection * (outerDirection.GetDotProduct(outerCentroid-outerIntercept));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::GetImpactParameters(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse)
{
    // sign convention for longitudinal distance:
    // -positive value means initial position is downstream of target position (i.e. cluster points back to vertex)
    transverse = initialDirection.GetCrossProduct(targetPosition-initialPosition).GetMagnitude();
    longitudinal = -initialDirection.GetDotProduct(targetPosition-initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArVertexHelper::m_numberOfLayersToFit = 20;
float LArVertexHelper::m_longitudinalWindow = 5.f; // cm

StatusCode LArVertexHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumberOfLayersToFit", m_numberOfLayersToFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LongitudinalWindow", m_longitudinalWindow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
