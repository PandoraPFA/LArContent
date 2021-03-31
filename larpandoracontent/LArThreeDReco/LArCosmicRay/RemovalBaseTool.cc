/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/RemovalBaseTool.cc
 *
 *  @brief  Implementation of the removal base class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h"

using namespace pandora;

namespace lar_content
{

RemovalBaseTool::RemovalBaseTool() :
    m_distanceToLine(0.5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool RemovalBaseTool::IsMuonEndpoint(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const bool ignoreHitType, const HitType &hitTypeToIgnore) const
{
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        if (ignoreHitType && (hitType == hitTypeToIgnore))
            continue;

        const Cluster *pMuonCluster(nullptr), *const pDeltaRayCluster(element.GetCluster(hitType));
    
        if (pAlgorithm->GetMuonCluster(element.GetOverlapResult().GetCommonMuonPfoList(), hitType, pMuonCluster) != STATUS_CODE_SUCCESS)
            return true;

        float xMinDR(-std::numeric_limits<float>::max()), xMaxDR(+std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanX(xMinDR, xMaxDR);
        
        float xMinCR(-std::numeric_limits<float>::max()), xMaxCR(+std::numeric_limits<float>::max());
        pMuonCluster->GetClusterSpanX(xMinCR, xMaxCR);

        if ((xMinDR < xMinCR) || (xMaxDR > xMaxCR))
            return true;
	
        float zMinDR(-std::numeric_limits<float>::max()), zMaxDR(+std::numeric_limits<float>::max());
        pDeltaRayCluster->GetClusterSpanZ(xMinDR, xMaxDR, zMinDR, zMaxDR);
        
        float zMinCR(-std::numeric_limits<float>::max()), zMaxCR(+std::numeric_limits<float>::max());
        pMuonCluster->GetClusterSpanZ(xMinCR, xMaxCR, zMinCR, zMaxCR);

        if ((zMinDR < zMinCR) || (zMaxDR > zMaxCR))
            return true;
    }

    return false;
}  
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool RemovalBaseTool::IsBestElement(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const HitType &hitType, 
    const TensorType::ElementList &elementList) const
{
    const float chiSquared(element.GetOverlapResult().GetReducedChi2());
    const unsigned int hitSum(element.GetClusterU()->GetNCaloHits() + element.GetClusterV()->GetNCaloHits() + element.GetClusterW()->GetNCaloHits());
            
    for (const TensorType::Element &testElement : elementList)
    {
        if (testElement.GetCluster(hitType) != element.GetCluster(hitType))
            continue;

        if ((testElement.GetClusterU() == element.GetClusterU()) && (testElement.GetClusterV() == element.GetClusterV()) && (testElement.GetClusterW() == element.GetClusterW()))
            continue;

        const unsigned int testHitSum(testElement.GetClusterU()->GetNCaloHits() + testElement.GetClusterV()->GetNCaloHits() + testElement.GetClusterW()->GetNCaloHits());
        const float testChiSquared(testElement.GetOverlapResult().GetReducedChi2());
        
        if ((testHitSum < hitSum) || ((testHitSum == hitSum) && (testChiSquared > chiSquared)))
            continue;

        if (this->PassElementChecks(pAlgorithm, testElement, hitType))
            return false;
    }

    return true;
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

        if (!this->IsCloseToLine(pCaloHit->GetPositionVector(), lowerBoundary, upperBoundary, m_distanceToLine))
            continue;

        collectedHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemovalBaseTool::ProjectDeltaRayPositions(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const Cluster *pCluster1(nullptr), *pCluster2(nullptr);
    
    for (const HitType &hitType1 : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if ((hitType2 == hitType) || (hitType1 == hitType2))
                continue;

            pCluster1 = element.GetCluster(hitType1);
            pCluster2 = element.GetCluster(hitType2);
            
            break;
        }
        break;
    }

    return pAlgorithm->GetProjectedPositions(pCluster1, pCluster2, projectedPositions);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode RemovalBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceToLine", m_distanceToLine));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
