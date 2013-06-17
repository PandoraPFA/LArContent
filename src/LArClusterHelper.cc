/**
 *  @file   LArClusterHelper.cc
 * 
 *  @brief  Implementation of the cluster helper class.
 * 
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"
#include "Helpers/XmlHelper.h"

#include "LArClusterHelper.h"

using namespace pandora;

namespace lar
{

float LArClusterHelper::GetClosestDistance(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2)
{
    const OrderedCaloHitList &orderedCaloHitList1(pCluster1->GetOrderedCaloHitList());
    const OrderedCaloHitList &orderedCaloHitList2(pCluster2->GetOrderedCaloHitList());

    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList1.begin(), iterEnd1 = orderedCaloHitList1.end(); iter1 != iterEnd1; ++iter1)
    {
        const CartesianVector position1(pCluster1->GetCentroid(iter1->first));
        
        for (OrderedCaloHitList::const_iterator iter2 = orderedCaloHitList2.begin(), iterEnd2 = orderedCaloHitList2.end(); iter2 != iterEnd2; ++iter2)
        {
            const CartesianVector position2(pCluster2->GetCentroid(iter2->first));
            const float distanceSquared((position2 - position1).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                distanceFound = true;
	    }
	}
    }

    if (distanceFound)
        return std::sqrt(closestDistanceSquared);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        const CartesianVector layerCentroid(pCluster->GetCentroid(iter->first));
        const float distanceSquared((layerCentroid - position).GetMagnitudeSquared());

        if (distanceSquared < closestDistanceSquared)
        {
            closestDistanceSquared = distanceSquared;
            distanceFound = true;
        }
    }
    
    if (distanceFound)
        return std::sqrt(closestDistanceSquared);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByInnerLayer(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int innerLayerLhs(pLhs->GetInnerPseudoLayer());
    const unsigned int innerLayerRhs(pRhs->GetInnerPseudoLayer());

    if( innerLayerLhs != innerLayerRhs )
      return (innerLayerLhs < innerLayerRhs);

    // Use SortByNOccupiedLayers method to resolve ties
    return SortByNOccupiedLayers(pLhs,pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNOccupiedLayers(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nOccupiedLayersLhs(pLhs->GetOrderedCaloHitList().size());
    const unsigned int nOccupiedLayersRhs(pRhs->GetOrderedCaloHitList().size());

    if (nOccupiedLayersLhs != nOccupiedLayersRhs)
        return (nOccupiedLayersLhs > nOccupiedLayersRhs);

    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNHits(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nHitsLhs(pLhs->GetNCaloHits());
    const unsigned int nHitsRhs(pRhs->GetNCaloHits());

    if (nHitsLhs != nHitsRhs)
        return (nHitsLhs > nHitsRhs);

    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsNode(const CartesianVector &parentVertex, const CartesianVector &daughterVertex)
{
    if ((parentVertex - daughterVertex).GetMagnitudeSquared() < m_maxNodeRadiusSquared)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsEmission(const LArPointingCluster::Vertex &parentPointingVertex, const pandora::CartesianVector &daughterVertex)
{
    return LArClusterHelper::IsPointing(daughterVertex, parentPointingVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsEmitted(const pandora::CartesianVector &parentVertex, const LArPointingCluster::Vertex &daughterPointingVertex)
{
    return LArClusterHelper::IsPointing(parentVertex, daughterPointingVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsPointing(const pandora::CartesianVector &vertex, const LArPointingCluster::Vertex &pointingVertex)
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

float LArClusterHelper::m_maxNodeRadiusSquared = 2.5f * 2.5f;
float LArClusterHelper::m_maxPointingLongitudinalDistance = 20.f;
float LArClusterHelper::m_minPointingLongitudinalDistance = -2.5f;
float LArClusterHelper::m_maxPointingTransverseDistance = 2.5f;
float LArClusterHelper::m_pointingAngularAllowance = 0.0175f; // tan (1 degree)

StatusCode LArClusterHelper::ReadSettings(const TiXmlHandle xmlHandle)
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
