/**<
 *  @file   larpandoracontent/LArHelpers/LArConnectionPathwayHelper.cc
 *
 *  @brief  Implementation of the connection pathway helper class.
 *
 *  $Log: $
 */

#include "Objects/CartesianVector.h"
#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"
#include "Pandora/Pandora.h"
#include "Pandora/Algorithm.h"

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"


#include <limits>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::SortByDistanceToPoint::operator()(const CartesianVector &lhs, const CartesianVector &rhs)
{
    return (m_referencePoint.GetDistanceSquared(lhs) < m_referencePoint.GetDistanceSquared(rhs));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::SortByDistanceToPoint::operator()(const CaloHit *const lhs, const CaloHit *const rhs)
{
    return (m_referencePoint.GetDistanceSquared(lhs->GetPositionVector()) < m_referencePoint.GetDistanceSquared(rhs->GetPositionVector()));
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Finding 3D Shower Start Position
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStarts3D(const pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerMatch &protoShowerMatch, const CartesianVector &nuVertexPosition, const float maxSeparationFromHit, 
    const float maxProjectionSeparation, const float maxXSeparation, CartesianPointVector &showerStarts3D)
{
    const ProtoShower &protoShowerU(protoShowerMatch.m_protoShowerU);
    const ProtoShower &protoShowerV(protoShowerMatch.m_protoShowerV);
    const ProtoShower &protoShowerW(protoShowerMatch.m_protoShowerW);
    const Consistency consistency(protoShowerMatch.m_consistencyType);

    bool uFound(false), vFound(false), wFound(false);
    CartesianVector uShowerStart3D(0.f, 0.f, 0.f), vShowerStart3D(0.f, 0.f, 0.f), wShowerStart3D(0.f, 0.f, 0.f);

    CaloHitList caloHitList3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHitList3D);

    if (consistency == Consistency::POSITION)
    {
        LArConnectionPathwayHelper::FindShowerStartFromPosition(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D);
        vShowerStart3D = uShowerStart3D;
        wShowerStart3D = uShowerStart3D;

        if (LArClusterHelper::GetClosestDistance(uShowerStart3D, caloHitList3D) < maxSeparationFromHit)
        {
            uFound = true; vFound = true; wFound = true;
        }
    }
    else if (consistency == Consistency::DIRECTION)
    {
        if (LArConnectionPathwayHelper::FindShowerStartFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, vShowerStart3D, wShowerStart3D))
        {
            if (LArClusterHelper::GetClosestDistance(uShowerStart3D, caloHitList3D) < maxSeparationFromHit)
                uFound = true;

            if (LArClusterHelper::GetClosestDistance(vShowerStart3D, caloHitList3D) < maxSeparationFromHit)
                vFound = true;

            if (LArClusterHelper::GetClosestDistance(wShowerStart3D, caloHitList3D) < maxSeparationFromHit)
                wFound = true;
        }
    }

    if (!uFound || (consistency == Consistency::X_PROJECTION))
    {
        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxProjectionSeparation, maxXSeparation, uShowerStart3D) || 
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxProjectionSeparation, maxXSeparation, uShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerU, protoShowerW, protoShowerV, maxProjectionSeparation, maxXSeparation, uShowerStart3D))
        {
                uFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxProjectionSeparation, maxXSeparation, vShowerStart3D) || 
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxProjectionSeparation, maxXSeparation, vShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerV, protoShowerU, protoShowerW, maxProjectionSeparation, maxXSeparation, vShowerStart3D))
        {
                vFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxProjectionSeparation, maxXSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxProjectionSeparation, maxXSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerW, protoShowerV, protoShowerU, maxProjectionSeparation, maxXSeparation, wShowerStart3D))
        {
                wFound = true;
        }
    }

    CartesianPointVector foundShowerStarts3D;

    if (uFound)
        foundShowerStarts3D.push_back(uShowerStart3D);

    if (vFound)
        foundShowerStarts3D.push_back(vShowerStart3D);

    if (wFound)
        foundShowerStarts3D.push_back(wShowerStart3D);

    std::sort(foundShowerStarts3D.begin(), foundShowerStarts3D.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertexPosition));

    if (foundShowerStarts3D.empty())
        return false;

    // ATTN: We need there to be three
    showerStarts3D.push_back(foundShowerStarts3D.front());
    showerStarts3D.push_back((foundShowerStarts3D.size() == 3) ? foundShowerStarts3D[1] : foundShowerStarts3D[0]);
    showerStarts3D.push_back(foundShowerStarts3D.back());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::FindShowerStartFromPosition(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    const CartesianVector showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector showerStartW(protoShowerW.m_showerCore.m_startPosition);

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV, showerStartW, showerStart3D, chi2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromDirection(const pandora::Algorithm *const pAlgorithm,
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &uShowerStart3D, 
    CartesianVector &vShowerStart3D, CartesianVector &wShowerStart3D)
{
    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    const CartesianVector &projectedNuVertexU(protoShowerU.m_connectionPathway.m_startPosition);
    const CartesianVector &projectedNuVertexV(protoShowerV.m_connectionPathway.m_startPosition);
    const CartesianVector &projectedNuVertexW(protoShowerW.m_connectionPathway.m_startPosition);

    const CartesianVector &peakDirectionU(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionV(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionW(protoShowerW.m_connectionPathway.m_startDirection);

    const CartesianVector &displacementU(showerStartU - projectedNuVertexU);
    const CartesianVector &displacementV(showerStartV - projectedNuVertexV);
    const CartesianVector &displacementW(showerStartW - projectedNuVertexW);

    const float transverseU(peakDirectionU.GetCrossProduct(displacementU).GetMagnitude());
    const float transverseV(peakDirectionV.GetCrossProduct(displacementV).GetMagnitude());
    const float transverseW(peakDirectionW.GetCrossProduct(displacementW).GetMagnitude());

    if ((transverseU > 1.f) || (transverseV > 1.f) || (transverseW > 1.f))
        return false;

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float cosThetaU(std::fabs(peakDirectionU.GetCosOpeningAngle(xAxis)));
    const float cosThetaV(std::fabs(peakDirectionV.GetCosOpeningAngle(xAxis)));
    const float cosThetaW(std::fabs(peakDirectionW.GetCosOpeningAngle(xAxis)));

    const float x1(showerStartU.GetX());
    const CartesianVector showerStartV1(projectedNuVertexV + (peakDirectionV * (std::fabs(x1 - projectedNuVertexV.GetX()) * cosThetaV)));
    const CartesianVector showerStartW1(projectedNuVertexW + (peakDirectionW * (std::fabs(x1 - projectedNuVertexW.GetX()) * cosThetaW)));

    const float x2(showerStartV.GetX());
    const CartesianVector showerStartU2(projectedNuVertexU + (peakDirectionU * (std::fabs(x2 - projectedNuVertexU.GetX()) * cosThetaU)));
    const CartesianVector showerStartW2(projectedNuVertexW + (peakDirectionW * (std::fabs(x2 - projectedNuVertexW.GetX()) * cosThetaW)));

    const float x3(showerStartW.GetX());
    const CartesianVector showerStartU3(projectedNuVertexU + (peakDirectionU * (std::fabs(x3 - projectedNuVertexU.GetX()) * cosThetaU)));
    const CartesianVector showerStartV3(projectedNuVertexV + (peakDirectionV * (std::fabs(x3 - projectedNuVertexV.GetX()) * cosThetaV)));

    float chi2(0.0);

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV1, showerStartW1, uShowerStart3D, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU2, showerStartV, showerStartW2, vShowerStart3D, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU3, showerStartV3, showerStartW, wShowerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromXProjection(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, const float maxXSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.m_spineHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_spineHitList.front()->GetHitType());
    const HitType hitType2(protoShower2.m_spineHitList.front()->GetHitType());

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower1, showerStart, maxXSeparation, showerStart1))
        return false;

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower2, showerStart, maxXSeparation, showerStart2))
        return false;

    float chi2(0.f);
    CartesianVector projection(0.f, 0.f, 0.f), projection1(0.f, 0.f, 0.f), projection2(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType1, hitType2, showerStart1, showerStart2, projection, chi2);
    const float separationSquared((projection - showerStart).GetMagnitudeSquared());

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType2, showerStart, showerStart2, projection1, chi2);
    const float separationSquared1((projection1 - showerStart1).GetMagnitudeSquared());

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);
    const float separationSquared2((projection2 - showerStart2).GetMagnitudeSquared());

    if ((separationSquared > maxSeparation * maxSeparation) || (separationSquared1 > maxSeparation * maxSeparation) || 
        (separationSquared2 > maxSeparation * maxSeparation))
    {
        return false;
    }

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), hitType, hitType1, hitType2, showerStart, showerStart1, showerStart2, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindClosestSpinePosition(const ProtoShower &protoShower, const CartesianVector &showerStart3D, const float maxXSeparation,
    CartesianVector &foundShowerStart)
{
    bool found(false);
    float lowestL(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
    {
        if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart3D.GetX()) > maxXSeparation)
            continue;

        float lVertex(protoShower.m_connectionPathway.m_startDirection.GetDotProduct(pCaloHit->GetPositionVector() - protoShower.m_connectionPathway.m_startPosition));

        if ((lVertex > 0.f) && (lVertex < lowestL))
        {
            lowestL = lVertex;
            foundShowerStart = pCaloHit->GetPositionVector();
            found = true;
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, const float maxXSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.m_spineHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_spineHitList.front()->GetHitType());

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower1, showerStart, maxXSeparation, showerStart1))
        return false;

    float chi2(0.f);
    CartesianVector projection2(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);

    // Now make sure that they agree...
    const CaloHitList &caloHitList2(protoShower2.m_spineHitList);
    const float separation(LArClusterHelper::GetClosestDistance(projection2, caloHitList2));

    if (separation > maxSeparation)
        return false;

    LArGeometryHelper::MergeTwoPositions3D(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Ordeing Functions
//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::GetMinMiddleMax(const float value1, const float value2, const float value3, float &minValue, float &middleValue, 
    float &maxValue)
{
    minValue = std::min(std::min(value1, value2), value3);
    maxValue = std::max(std::max(value1, value2), value3);
    middleValue = minValue;

    for (const float value : {value1, value2, value3})
    {
        if ((std::fabs(value - minValue) < std::numeric_limits<float>::epsilon()) || 
            (std::fabs(value - maxValue) < std::numeric_limits<float>::epsilon()))
        {
            continue;
        }

        middleValue = value;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
