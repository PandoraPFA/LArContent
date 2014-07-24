/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

void BoundedClusterMergingAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters,
    const ClusterToListNameMap &clusterToListNameMap) const
{
    ClusterAssociationMap clusterAssociationMap;

    for (ClusterList::const_iterator pIter = pfoClusters.begin(), pIterEnd = pfoClusters.end(); pIter != pIterEnd; ++pIter)
    {
        Cluster *pPfoCluster(*pIter);
        const TwoDSlidingShowerFitResult fitResult(pPfoCluster, m_slidingFitWindow);

        ShowerPositionMap showerPositionMap;
        const XSampling xSampling(fitResult.GetShowerFitResult());
        this->GetShowerPositionMap(fitResult, xSampling, showerPositionMap);

        for (ClusterList::const_iterator rIter = remnantClusters.begin(), rIterEnd = remnantClusters.end(); rIter != rIterEnd; ++rIter)
        {
            Cluster *pRemnantCluster(*rIter);
            const float boundedFraction(this->GetBoundedFraction(pRemnantCluster, xSampling, showerPositionMap));

            if (boundedFraction < m_minBoundedFraction)
                continue;

            AssociationDetails &associationDetails(clusterAssociationMap[pRemnantCluster]);

            if (!associationDetails.insert(AssociationDetails::value_type(pPfoCluster, boundedFraction)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }

    this->MakeClusterMerges(clusterAssociationMap, clusterToListNameMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::GetShowerPositionMap(const TwoDSlidingShowerFitResult &fitResult, const XSampling &xSampling,
    ShowerPositionMap &showerPositionMap) const
{
    for (float x = xSampling.m_minX; x < xSampling.m_maxX; x += xSampling.m_xPitch)
    {
        FloatVector edgePositions;
        fitResult.GetShowerEdges(x, edgePositions);

        if (edgePositions.size() < 2)
            continue;

        std::sort(edgePositions.begin(), edgePositions.end());
        const int xBin((x - xSampling.m_minX) / xSampling.m_xPitch);
        showerPositionMap.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, edgePositions.front(), edgePositions.back())));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float BoundedClusterMergingAlgorithm::GetBoundedFraction(const Cluster *const pCluster, const XSampling &xSampling, const ShowerPositionMap &showerPositionMap) const
{
    if (((xSampling.m_maxX - xSampling.m_minX) < std::numeric_limits<float>::epsilon()) ||
        (xSampling.m_xPitch < std::numeric_limits<float>::epsilon()) || (0 == pCluster->GetNCaloHits()))
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    unsigned int nMatchedHits(0);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            if ((x < xSampling.m_minX) || (x > xSampling.m_maxX))
                continue;

            const int xBin((x - xSampling.m_minX) / xSampling.m_xPitch);

            ShowerPositionMap::const_iterator positionIter = showerPositionMap.find(xBin);

            if ((showerPositionMap.end() != positionIter) && (z > positionIter->second.GetLowEdgeZ()) && (z < positionIter->second.GetHighEdgeZ()))
                ++nMatchedHits;
        }
    }

    return (static_cast<float>(nMatchedHits) / static_cast<float>(pCluster->GetNCaloHits()));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BoundedClusterMergingAlgorithm::XSampling::XSampling(const TwoDSlidingFitResult &fitResult)
{
    fitResult.GetMinAndMaxX(m_minX, m_maxX);
    const int nPoints(fitResult.GetMaxLayer() - fitResult.GetMinLayer());

    if (0 >= nPoints)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    m_xPitch = (m_maxX - m_minX) / static_cast<float>(nPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BoundedClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "SlidingFitWindow", m_slidingFitWindow));

    m_minBoundedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinBoundedFraction", m_minBoundedFraction));

    return ClusterMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
