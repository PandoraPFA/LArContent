/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/RemovalBaseTool.cc
 *
 *  @brief  Implementation of the removal base class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

RemovalBaseTool::RemovalBaseTool() :
    m_xOverlapWindow(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemovalBaseTool::GetMuonCluster(const TensorType::Element &element, const HitType &hitType, const Cluster *&pMuonCluster) const
{
    const PfoList commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());

    if (commonMuonPfoList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(commonMuonPfoList.front(), hitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    pMuonCluster = muonClusterList.front();

    return STATUS_CODE_SUCCESS;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

bool RemovalBaseTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, const CartesianVector &lineEnd, const float distanceToLine) const
{
    CartesianVector lineDirection(lineStart - lineEnd);
    lineDirection = lineDirection.GetUnitVector();
    
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
       return false;

    return true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

bool RemovalBaseTool::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float segmentBoundaryGradient = (-1.f) * (upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ()) / segmentBoundaryGradient + upperBoundary.GetX());
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ()) / segmentBoundaryGradient + lowerBoundary.GetX());
    
    if (std::fabs(xPointOnUpperLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;

    if (std::fabs(xPointOnLowerLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;
    
    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RemovalBaseTool::FindExtrapolatedHits(const Cluster *const pCluster, const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary,
    CaloHitList &collectedHits) const
{
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (!this->IsInLineSegment(lowerBoundary, upperBoundary, pCaloHit->GetPositionVector()))
            continue;

        if (!this->IsCloseToLine(pCaloHit->GetPositionVector(), lowerBoundary, upperBoundary, 0.5))
            continue;

        collectedHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

float RemovalBaseTool::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CaloHit *const pTestCaloHit : caloHitList)
    {
        const CartesianVector &position(pTestCaloHit->GetPositionVector());
        float separationSquared((position - referencePoint).GetMagnitudeSquared());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return std::sqrt(shortestDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float RemovalBaseTool::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        const float separationSquared((testPosition - referencePoint).GetMagnitudeSquared());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return std::sqrt(shortestDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector RemovalBaseTool::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector, const Cluster *const pCluster) const
{
    CartesianVector closestPoint(0.f,0.f,0.f);
    float shortestDistanceSquared(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > 0.5f)
            continue;

        const float separationSquared((testPosition - referencePoint).GetMagnitude());

        if (separationSquared > 5.f)
            continue;

        if (separationSquared < shortestDistanceSquared)
        {
            shortestDistanceSquared = separationSquared;
            closestPoint = testPosition;
        }
    }

    return closestPoint;
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void RemovalBaseTool::ProjectMuonPositions(const TensorType::Element &element, const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            const Cluster *pMuonCluster1(nullptr), *pMuonCluster2(nullptr);
            
            if (this->GetMuonCluster(element, hitType1, pMuonCluster1) != STATUS_CODE_SUCCESS)
                return;
                        
            if (this->GetMuonCluster(element, hitType2, pMuonCluster2) != STATUS_CODE_SUCCESS)
                return;

            this->ProjectPositions(pMuonCluster1, pMuonCluster2, projectedPositions);

            break;
        }
        
        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RemovalBaseTool::ProjectDeltaRayPositions(const TensorType::Element &element, const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            this->ProjectPositions(element.GetCluster(hitType1), element.GetCluster(hitType2), projectedPositions);

            break;
        }
        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RemovalBaseTool::ProjectPositions(const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianPointVector &projectedPositions) const
{    
    float xMin1(-std::numeric_limits<float>::max()), xMax1(+std::numeric_limits<float>::max());
    float xMin2(-std::numeric_limits<float>::max()), xMax2(+std::numeric_limits<float>::max());

    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, xMin2) - xPitch);
    const float xMax(std::min(xMax1, xMax2) + xPitch);
    const float xOverlap(xMax - xMin);

    // this has already been done...
    if (xOverlap < std::numeric_limits<float>::epsilon())
      return;
    
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    if (hitType1 == hitType2)
        return;

    const unsigned int nPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    // Projection into third view
    for (unsigned int n = 0; n < nPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMin1(0.f), zMin2(0.f), zMax1(0.f), zMax2(0.f);
            pCluster1->GetClusterSpanZ(xmin, xmax, zMin1, zMax1);
            pCluster2->GetClusterSpanZ(xmin, xmax, zMin2, zMax2);

            const float z1(0.5f * (zMin1 + zMax1));
            const float z2(0.5f * (zMin2 + zMax2));

            // could make use of the chi-squared?
            float chi2(0.f);
            CartesianVector projection(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, CartesianVector(x, 0.f, z1), CartesianVector(x, 0.f, z2), projection, chi2);

            projectedPositions.push_back(projection);
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (statusCodeException.GetStatusCode() != STATUS_CODE_NOT_FOUND)
                throw statusCodeException.GetStatusCode();

            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemovalBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapWindow", m_xOverlapWindow));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
