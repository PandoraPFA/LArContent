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
// Consistency checks
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation)
{
    float metric(0.0);

    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    const float xSeparationUV(std::fabs(showerStartU.GetX() - showerStartV.GetX()));
    const float xSeparationUW(std::fabs(showerStartU.GetX() - showerStartW.GetX()));
    const float xSeparationVW(std::fabs(showerStartV.GetX() - showerStartW.GetX()));

    if ((xSeparationUV > maxXSeparation) || (xSeparationUW > maxXSeparation) || (xSeparationVW > maxXSeparation))
        return false;

    float chi2(0.f);
    CartesianVector projectionU(0.f, 0.f, 0.f), projectionV(0.f, 0.f, 0.f), projectionW(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, showerStartV, showerStartW, projectionU, chi2);
    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, showerStartW, showerStartU, projectionV, chi2);
    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, showerStartU, showerStartV, projectionW, chi2);

    const float separationU((projectionU - showerStartU).GetMagnitude());
    const float separationV((projectionV - showerStartV).GetMagnitude());
    const float separationW((projectionW - showerStartW).GetMagnitude());

    metric = (separationU + separationV + separationW) / 3.f;

    return (metric < maxSeparation);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle)
{
    const CartesianVector &directionU1(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &directionV1(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &directionW1(protoShowerW.m_connectionPathway.m_startDirection);

    if (LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU1, directionV1, directionW1, maxOpeningAngle))
    {
        return true;
    }
    else
    {
        const bool isDownstream(protoShowerW.m_showerCore.m_startPosition.GetZ() > protoShowerW.m_connectionPathway.m_startPosition.GetZ());

        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
            spinePositionsU.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
            spinePositionsV.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
            spinePositionsW.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const CartesianVector directionU2(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionV2(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionW2(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        return LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU2, directionV2, directionW2, maxOpeningAngle);
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, CartesianVector directionU, 
    CartesianVector directionV, CartesianVector directionW, const float maxOpeningAngle)
{
    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    float radians((2.f * M_PI) / 180.f);
    bool isIsochronous((wireDeviationU < radians) && (wireDeviationV < radians) && (wireDeviationW < radians));

    if (isIsochronous)
    {
        int positiveCount(0);
        positiveCount += directionU.GetX() > 0.f ? 0 : 1;
        positiveCount += directionV.GetX() > 0.f ? 0 : 1;
        positiveCount += directionW.GetX() > 0.f ? 0 : 1;

        if (positiveCount >= 2)
        {
            directionU = CartesianVector(std::fabs(directionU.GetX()), 0.f, directionU.GetZ());
            directionV = CartesianVector(std::fabs(directionV.GetX()), 0.f, directionV.GetZ());
            directionW = CartesianVector(std::fabs(directionW.GetX()), 0.f, directionW.GetZ());
        }
    }

    if (directionU.GetX() * directionV.GetX() < 0.f)
        return false;

    if (directionU.GetX() * directionW.GetX() < 0.f)
        return false;

    if (directionV.GetX() * directionW.GetX() < 0.f)
        return false;
    
    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    if (isIsochronous)
    {
        openingAngleU = std::min(openingAngleU, 180.f - openingAngleU);
        openingAngleV = std::min(openingAngleV, 180.f - openingAngleV);
        openingAngleW = std::min(openingAngleW, 180.f - openingAngleW);
    }

    if ((openingAngleU > maxOpeningAngle) || (openingAngleV > maxOpeningAngle) || (openingAngleW > maxOpeningAngle))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Finding 3D Shower Start Position
//------------------------------------------------------------------------------------------------------------------------------------------
// maxSeparationFromHit = 3.f, maxProjectionSeparation = 5.f
bool LArConnectionPathwayHelper::FindShowerStarts3D(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerMatch &protoShowerMatch, const CartesianVector &nuVertexPosition, const float maxSeparationFromHit, 
    const float maxProjectionSeparation, CartesianPointVector &showerStarts3D)
{
    const ElectronProtoShower &protoShowerU(protoShowerMatch.m_protoShowerU);
    const ElectronProtoShower &protoShowerV(protoShowerMatch.m_protoShowerV);
    const ElectronProtoShower &protoShowerW(protoShowerMatch.m_protoShowerW);
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
        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxProjectionSeparation, uShowerStart3D) || 
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxProjectionSeparation, uShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerU, protoShowerW, protoShowerV, maxProjectionSeparation, uShowerStart3D))
        {
                uFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxProjectionSeparation, vShowerStart3D) || 
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxProjectionSeparation, vShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerV, protoShowerU, protoShowerW, maxProjectionSeparation, vShowerStart3D))
        {
                vFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerStartFromXProjection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxProjectionSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxProjectionSeparation, wShowerStart3D) ||
            LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pAlgorithm, protoShowerW, protoShowerV, protoShowerU, maxProjectionSeparation, wShowerStart3D))
        {
                wFound = true;
        }
    }

    CartesianPointVector tempShowerStarts3D;

    if (uFound)
        tempShowerStarts3D.push_back(uShowerStart3D);

    if (vFound)
        tempShowerStarts3D.push_back(vShowerStart3D);

    if (wFound)
        tempShowerStarts3D.push_back(wShowerStart3D);

    std::sort(tempShowerStarts3D.begin(), tempShowerStarts3D.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertexPosition));

    if (tempShowerStarts3D.empty())
        return false;

    showerStarts3D.push_back(tempShowerStarts3D.front());
    showerStarts3D.push_back((tempShowerStarts3D.size() == 3) ? tempShowerStarts3D[1] : tempShowerStarts3D[0]);
    showerStarts3D.push_back(tempShowerStarts3D.back());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromPosition(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    const CartesianVector showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector showerStartW(protoShowerW.m_showerCore.m_startPosition);

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV, showerStartW, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromDirection(pandora::Algorithm *const pAlgorithm, const CartesianVector nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    CartesianVector uShowerStart3D(0.f, 0.f, 0.f), vShowerStart3D(0.f, 0.f, 0.f), wShowerStart3D(0.f, 0.f, 0.f);

    if(!LArConnectionPathwayHelper::FindShowerStartFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, 
        vShowerStart3D, wShowerStart3D))
    {
        return false;
    }

    float minSeparation(std::numeric_limits<float>::max());

    for (const CartesianVector showerStart : {uShowerStart3D, vShowerStart3D, wShowerStart3D})
    {
        const float separation((showerStart - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            minSeparation = separation;
            showerStart3D = showerStart;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStartFromDirection(pandora::Algorithm *const pAlgorithm,
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

bool LArConnectionPathwayHelper::FindShowerStartFromXProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.m_spineHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_spineHitList.front()->GetHitType());
    const HitType hitType2(protoShower2.m_spineHitList.front()->GetHitType());

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower1, showerStart, showerStart1))
        return false;

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower2, showerStart, showerStart2))
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

bool LArConnectionPathwayHelper::FindClosestSpinePosition(const ProtoShower &protoShower, const CartesianVector &showerStart3D, 
    CartesianVector &foundShowerStart)
{
    bool found(false);
    float lowestL(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
    {
        if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart3D.GetX()) > 0.5f)
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

bool LArConnectionPathwayHelper::FindShowerStartFromXProjectionRelaxed(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.m_spineHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_spineHitList.front()->GetHitType());

    if (!LArConnectionPathwayHelper::FindClosestSpinePosition(protoShower1, showerStart, showerStart1))
        return false;

    // Now make sure that they agree...
    const CaloHitList &caloHitList2(protoShower2.m_spineHitList);

    float chi2(0.f);
    CartesianVector projection2(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);
    const float separation(LArClusterHelper::GetClosestDistance(projection2, caloHitList2));

    if (separation > maxSeparation)
        return false;

    LArGeometryHelper::MergeTwoPositions3D(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, showerStart3D, chi2);

    return true;
}

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


//------------------------------------------------------------------------------------------------------------------------------------------
// Sorting Functions
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

} // namespace lar_content
