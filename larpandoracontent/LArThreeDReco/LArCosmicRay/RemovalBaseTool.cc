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

RemovalBaseTool::RemovalBaseTool()
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

bool RemovalBaseTool::IsBestElement(const TensorType::Element &element, const HitType &hitType, const TensorType::ElementList &elementList) const
{
    float chiSquared(element.GetOverlapResult().GetReducedChi2());
    unsigned int hitNumber(element.GetClusterU()->GetNCaloHits() + element.GetClusterV()->GetNCaloHits() + element.GetClusterW()->GetNCaloHits());
            
    for (const TensorType::Element &testElement : elementList)
    {
        if (testElement.GetCluster(hitType) != element.GetCluster(hitType))
            continue;

        if ((testElement.GetClusterU() == element.GetClusterU()) && (testElement.GetClusterV() == element.GetClusterV()) && (testElement.GetClusterW() == element.GetClusterW()))
            continue;

        const unsigned int hitSum(testElement.GetClusterU()->GetNCaloHits() + testElement.GetClusterV()->GetNCaloHits() + testElement.GetClusterW()->GetNCaloHits());

        if ((hitSum == hitNumber) && (testElement.GetOverlapResult().GetReducedChi2() < chiSquared))
            return false;
        
        if (hitSum > hitNumber)
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

        if (!this->IsCloseToLine(pCaloHit->GetPositionVector(), lowerBoundary, upperBoundary, 0.5))
            continue;

        collectedHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemovalBaseTool::ProjectDeltaRayPositions(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const HitType &hitType, CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    const Cluster *pCluster1(nullptr), *pCluster2(nullptr);
    
    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == hitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
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

StatusCode RemovalBaseTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
