/**
 *  @file   larpandoracontent/LArHelpers/LArConnectionPathwayHelper.cc
 *
 *  @brief  Implementation of the connection pathway helper class.
 *
 *  $Log: $
 */

#include "Pandora/Algorithm.h"
#include "Pandora/Pandora.h"
#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

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

bool LArConnectionPathwayHelper::FindShowerStarts3D(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const ProtoShowerMatch &protoShowerMatch, const CartesianVector &nuVertexPosition, const float maxSeparationFromHit,
    const float maxProjectionSeparation, const float maxXSeparation, CartesianPointVector &showerStarts3D)
{
    const ProtoShower &protoShowerU(protoShowerMatch.GetProtoShowerU());
    const ProtoShower &protoShowerV(protoShowerMatch.GetProtoShowerV());
    const ProtoShower &protoShowerW(protoShowerMatch.GetProtoShowerW());
    const Consistency consistency(protoShowerMatch.GetConsistencyType());

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
            uFound = true;
            vFound = true;
            wFound = true;
        }
    }
    else if (consistency == Consistency::DIRECTION)
    {
        if (LArConnectionPathwayHelper::FindShowerStartFromDirection(
                pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, vShowerStart3D, wShowerStart3D))
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
        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(
                pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxProjectionSeparation, maxXSeparation, uShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(
                pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxProjectionSeparation, maxXSeparation, uShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(
                pAlgorithm, protoShowerU, protoShowerW, protoShowerV, maxProjectionSeparation, maxXSeparation, uShowerStart3D))
        {
            uFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(
                pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxProjectionSeparation, maxXSeparation, vShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(
                pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxProjectionSeparation, maxXSeparation, vShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(
                pAlgorithm, protoShowerV, protoShowerU, protoShowerW, maxProjectionSeparation, maxXSeparation, vShowerStart3D))
        {
            vFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(
                pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxProjectionSeparation, maxXSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(
                pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxProjectionSeparation, maxXSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(
                pAlgorithm, protoShowerW, protoShowerV, protoShowerU, maxProjectionSeparation, maxXSeparation, wShowerStart3D))
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

    if (foundShowerStarts3D.empty())
        return false;

    std::sort(foundShowerStarts3D.begin(), foundShowerStarts3D.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertexPosition));

    // ATTN: We need there to be three
    showerStarts3D.push_back(foundShowerStarts3D.front());
    showerStarts3D.push_back((foundShowerStarts3D.size() == 3) ? foundShowerStarts3D[1] : foundShowerStarts3D[0]);
    showerStarts3D.push_back(foundShowerStarts3D.back());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::FindShowerStartFromPosition(const Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    const CartesianVector showerStartU(protoShowerU.GetShowerCore().GetStartPosition());
    const CartesianVector showerStartV(protoShowerV.GetShowerCore().GetStartPosition());
    const CartesianVector showerStartW(protoShowerW.GetShowerCore().GetStartPosition());

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(
        pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV, showerStartW, showerStart3D, chi2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromDirection(const Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &uShowerStart3D, CartesianVector &vShowerStart3D,
    CartesianVector &wShowerStart3D)
{

    if (!LArConnectionPathwayHelper::FindShowerStartFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D))
        return false;

    if (!LArConnectionPathwayHelper::FindShowerStartFromDirection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, vShowerStart3D))
        return false;

    if (!LArConnectionPathwayHelper::FindShowerStartFromDirection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, wShowerStart3D))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromDirection(const Algorithm *const pAlgorithm, const ProtoShower &protoShower,
    const ProtoShower &protoShowerA, const ProtoShower &protoShowerB, CartesianVector &showerStart3D)
{
    const CartesianVector &showerStart(protoShower.GetShowerCore().GetStartPosition());
    const float x(showerStart.GetX());
    CartesianVector showerStartA(0.f, 0.f, 0.f), showerStartB(0.f, 0.f, 0.f);

    if (!LArConnectionPathwayHelper::ProjectShowerStartByDirection(protoShowerA, x, showerStartA))
        return false;

    if (!LArConnectionPathwayHelper::ProjectShowerStartByDirection(protoShowerB, x, showerStartB))
        return false;

    float chi2(0.0);
    const HitType hitType(protoShower.GetSpineHitList().front()->GetHitType());
    const HitType hitTypeA(protoShowerA.GetSpineHitList().front()->GetHitType());
    const HitType hitTypeB(protoShowerB.GetSpineHitList().front()->GetHitType());

    LArGeometryHelper::MergeThreePositions3D(
        pAlgorithm->GetPandora(), hitType, hitTypeA, hitTypeB, showerStart, showerStartA, showerStartB, showerStart3D, chi2);

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::ProjectShowerStartByDirection(const ProtoShower &protoShower, const float x, CartesianVector &showerStart2D)
{
    const CartesianVector &viewShowerStart(protoShower.GetShowerCore().GetStartPosition());
    const CartesianVector &viewNuVertex(protoShower.GetConnectionPathway().GetStartPosition());
    const CartesianVector &viewPeakDirection(protoShower.GetConnectionPathway().GetStartDirection());
    const CartesianVector &displacement(viewShowerStart - viewNuVertex);
    const float transverse(viewPeakDirection.GetCrossProduct(displacement).GetMagnitude());

    if (transverse > 1.f)
        return false;

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float cosTheta(std::fabs(viewPeakDirection.GetCosOpeningAngle(xAxis)));

    showerStart2D = viewNuVertex + (viewPeakDirection * (std::fabs(x - viewNuVertex.GetX()) * cosTheta));

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromXProjection(const Algorithm *const pAlgorithm, const ProtoShower &protoShower,
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, const float maxXSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.GetShowerCore().GetStartPosition());
    CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.GetSpineHitList().front()->GetHitType());
    const HitType hitType1(protoShower1.GetSpineHitList().front()->GetHitType());
    const HitType hitType2(protoShower2.GetSpineHitList().front()->GetHitType());

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

    LArGeometryHelper::MergeThreePositions3D(
        pAlgorithm->GetPandora(), hitType, hitType1, hitType2, showerStart, showerStart1, showerStart2, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindClosestSpinePosition(
    const ProtoShower &protoShower, const CartesianVector &showerStart3D, const float maxXSeparation, CartesianVector &foundShowerStart)
{
    bool found(false);
    float lowestL(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower.GetSpineHitList())
    {
        if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart3D.GetX()) > maxXSeparation)
            continue;

        const float lVertex(protoShower.GetConnectionPathway().GetStartDirection().GetDotProduct(
            pCaloHit->GetPositionVector() - protoShower.GetConnectionPathway().GetStartPosition()));

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

bool LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(const Algorithm *const pAlgorithm,
    const ProtoShower &protoShower, const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation,
    const float maxXSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.GetShowerCore().GetStartPosition());
    CartesianVector showerStart1(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.GetSpineHitList().front()->GetHitType());
    const HitType hitType1(protoShower1.GetSpineHitList().front()->GetHitType());

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower1, showerStart, maxXSeparation, showerStart1))
        return false;

    float chi2(0.f);
    CartesianVector projection2(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);

    // Now make sure that they agree...
    const CaloHitList &caloHitList2(protoShower2.GetSpineHitList());
    const float separation(LArClusterHelper::GetClosestDistance(projection2, caloHitList2));

    if (separation > maxSeparation)
        return false;

    LArGeometryHelper::MergeTwoPositions3D(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Ordeing Functions
//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::GetMinMiddleMax(
    const float value1, const float value2, const float value3, float &minValue, float &middleValue, float &maxValue)
{
    FloatVector values({value1, value2, value3});

    std::sort(values.begin(), values.end());

    minValue = values.at(0);
    middleValue = values.at(1);
    maxValue = values.at(2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
