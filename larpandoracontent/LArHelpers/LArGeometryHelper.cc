/**
 *  @file   larpandoracontent/LArHelpers/LArGeometryHelper.cc
 *
 *  @brief  Implementation of the geometry helper class.
 *
 *  $Log: $
 */

#include "Managers/PluginManager.h"
#include "Managers/GeometryManager.h"

#include "Objects/CartesianVector.h"
#include "Objects/DetectorGap.h"

#include "Pandora/Pandora.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

using namespace pandora;

namespace lar_content
{

LArGeometryHelper::PseudoLayerInstanceMap LArGeometryHelper::m_pseudolayerInstanceMap;
LArGeometryHelper::TransformationInstanceMap LArGeometryHelper::m_transformationInstanceMap;

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::MergeTwoPositions(const Pandora &pandora, const HitType view1, const HitType view2, const float position1, const float position2)
{
    if (view1 == view2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        return LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoW(position1, position2);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        return LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoW(position2, position1);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        return LArGeometryHelper::GetLArTransformationPlugin(pandora)->WUtoV(position1, position2);
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        return LArGeometryHelper::GetLArTransformationPlugin(pandora)->WUtoV(position2, position1);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        return LArGeometryHelper::GetLArTransformationPlugin(pandora)->VWtoU(position1, position2);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        return LArGeometryHelper::GetLArTransformationPlugin(pandora)->VWtoU(position2, position1);
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::MergeTwoDirections(const Pandora &pandora, const HitType view1, const HitType view2,
    const CartesianVector &direction1, const CartesianVector &direction2)
{
    // x components must be consistent
    if (direction1.GetX() * direction2.GetX() < 0.f)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // scale x components to a common value
    const float s1((std::fabs(direction1.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction2.GetX()) : 1.f);
    const float s2((std::fabs(direction2.GetX()) > std::numeric_limits<float>::epsilon()) ? 100.f * std::fabs(direction1.GetX()) : 1.f);

    float pX(s1 * direction1.GetX()), pU(0.f), pV(0.f), pW(0.f);
    HitType newView(HIT_CUSTOM);

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        pU = s1 * direction1.GetZ(); pV = s2 * direction2.GetZ(); newView = TPC_VIEW_W;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        pV = s1 * direction1.GetZ(); pU = s2 * direction2.GetZ(); newView = TPC_VIEW_W;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        pW = s1 * direction1.GetZ(); pU = s2 * direction2.GetZ(); newView = TPC_VIEW_V;
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        pU = s1 * direction1.GetZ(); pW = s2 * direction2.GetZ(); newView = TPC_VIEW_V;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        pV = s1 * direction1.GetZ(); pW = s2 * direction2.GetZ(); newView = TPC_VIEW_U;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        pW = s1 * direction1.GetZ(); pV = s2 * direction2.GetZ(); newView = TPC_VIEW_U;
    }

    if (newView == TPC_VIEW_W)
        return CartesianVector(pX, 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->PUPVtoPW(pU, pV)).GetUnitVector();

    if (newView == TPC_VIEW_U)
        return CartesianVector(pX, 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->PVPWtoPU(pV, pW)).GetUnitVector();

    if (newView == TPC_VIEW_V)
        return CartesianVector(pX, 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->PWPUtoPV(pW, pU)).GetUnitVector();

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions(const Pandora &pandora, const HitType view1, const HitType view2, const CartesianVector &position1,
    const CartesianVector &position2, CartesianVector &position3, float& chiSquared)
{
    if (view1 == view2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float X3((position1.GetX() + position2.GetX() ) / 2.f);
    const float Z1(position1.GetZ());
    const float Z2(position2.GetZ());
    const float Z3(LArGeometryHelper::MergeTwoPositions(pandora, view1, view2, Z1, Z2));

    position3.SetValues(X3, 0.f, Z3);
    const float sigmaUVW(LArGeometryHelper::GetLArTransformationPlugin(pandora)->GetSigmaUVW());
    chiSquared = ((X3 - position1.GetX()) * (X3 - position1.GetX()) + (X3 - position2.GetX()) * (X3 - position2.GetX())) / (sigmaUVW * sigmaUVW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions(const Pandora &pandora, const HitType view1, const HitType view2, const CartesianVector &position1,
    const CartesianVector &position2, CartesianVector &outputU, CartesianVector &outputV, CartesianVector &outputW, float& chiSquared)
{
    if (view1 == view2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector position3(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(pandora, view1, view2, position1, position2, position3, chiSquared);

    float aveU(0.f), aveV(0.f), aveW(0.f);
    const float Z1(position1.GetZ()), Z2(position2.GetZ()), Z3(position3.GetZ()), aveX(position3.GetX());

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        aveU = Z1; aveV = Z2; aveW = Z3;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        aveV = Z1; aveW = Z2; aveU = Z3;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        aveW = Z1; aveU = Z2; aveV = Z3;
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        aveV = Z1; aveU = Z2; aveW = Z3;
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        aveW = Z1; aveV = Z2; aveU = Z3;
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        aveU = Z1; aveW = Z2; aveV = Z3;
    }

    outputU.SetValues( aveX, 0.f, aveU );
    outputV.SetValues( aveX, 0.f, aveV );
    outputW.SetValues( aveX, 0.f, aveW );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions(const Pandora &pandora, const HitType view1, const HitType view2, const HitType view3,
    const CartesianVector &position1, const CartesianVector &position2, const CartesianVector &position3, CartesianVector &outputU,
    CartesianVector &outputV, CartesianVector &outputW, float &chiSquared)
{
    if ((view1 == view2) || (view2 == view3) || (view3 == view1))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_V))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position1, position2, position3, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_W))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position3, position1, position2, outputU, outputV, outputW, chiSquared );
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_U))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position2, position3, position1, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_V) && (view2 == TPC_VIEW_U))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position2, position1, position3, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_W) && (view2 == TPC_VIEW_V))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position3, position2, position1, outputU, outputV, outputW, chiSquared);
    }

    if ((view1 == TPC_VIEW_U) && (view2 == TPC_VIEW_W))
    {
        return LArGeometryHelper::MergeThreePositions(pandora, position1, position3, position2, outputU, outputV, outputW, chiSquared);
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions(const Pandora &pandora, const CartesianVector &positionU, const CartesianVector &positionV,
    const CartesianVector &positionW, CartesianVector &outputU, CartesianVector &outputV, CartesianVector &outputW, float &chiSquared)
{
    const float YfromUV(LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoY(positionU.GetZ(), positionV.GetZ()));
    const float ZfromUV(LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoZ(positionU.GetZ(), positionV.GetZ()));

    const float aveX((positionU.GetX() + positionV.GetX() + positionW.GetX()) / 3.f);
    const float aveY(YfromUV);
    const float aveW((positionW.GetZ() + 2.f * ZfromUV ) / 3.f);

    const float aveU(LArGeometryHelper::GetLArTransformationPlugin(pandora)->YZtoU(aveY, aveW));
    const float aveV(LArGeometryHelper::GetLArTransformationPlugin(pandora)->YZtoV(aveY, aveW));

    outputU.SetValues(aveX, 0.f, aveU);
    outputV.SetValues(aveX, 0.f, aveV);
    outputW.SetValues(aveX, 0.f, aveW);

    const float sigmaUVW(LArGeometryHelper::GetLArTransformationPlugin(pandora)->GetSigmaUVW());
    chiSquared = ((outputU.GetX() - positionU.GetX()) * (outputU.GetX() - positionU.GetX()) +
        (outputV.GetX() - positionV.GetX()) * (outputV.GetX() - positionV.GetX()) +
        (outputW.GetX() - positionW.GetX()) * (outputW.GetX() - positionW.GetX()) +
        (outputU.GetZ() - positionU.GetZ()) * (outputU.GetZ() - positionU.GetZ()) +
        (outputV.GetZ() - positionV.GetZ()) * (outputV.GetZ() - positionV.GetZ()) +
        (outputW.GetZ() - positionW.GetZ()) * (outputW.GetZ() - positionW.GetZ())) / (sigmaUVW * sigmaUVW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeTwoPositions3D(const Pandora &pandora, const HitType view1, const HitType view2, const CartesianVector &position1,
    const CartesianVector &position2, CartesianVector &position3D, float &chiSquared)
{
    CartesianVector positionU(0.f, 0.f, 0.f), positionV(0.f, 0.f, 0.f), positionW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(pandora, view1, view2, position1, position2, positionU, positionV, positionW, chiSquared);

    position3D.SetValues(positionW.GetX(), LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoY(positionU.GetZ(),positionV.GetZ()),
      LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoZ(positionU.GetZ(),positionV.GetZ()) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGeometryHelper::MergeThreePositions3D(const Pandora &pandora, const HitType view1, const HitType view2, const HitType view3,
    const CartesianVector &position1, const CartesianVector &position2, const CartesianVector &position3, CartesianVector &position3D, float &chiSquared)
{
    CartesianVector positionU(0.f, 0.f, 0.f), positionV(0.f, 0.f, 0.f), positionW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeThreePositions(pandora, view1, view2, view3, position1, position2, position3, positionU, positionV, positionW, chiSquared);

    position3D.SetValues(positionW.GetX(), LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoY(positionU.GetZ(),positionV.GetZ()),
        LArGeometryHelper::GetLArTransformationPlugin(pandora)->UVtoZ(positionU.GetZ(),positionV.GetZ()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::ProjectPosition(const Pandora &pandora, const CartesianVector &position3D, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        return CartesianVector(position3D.GetX(), 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->YZtoU(position3D.GetY(), position3D.GetZ()));
    }

    else if (view == TPC_VIEW_V)
    {
        return CartesianVector(position3D.GetX(), 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->YZtoV(position3D.GetY(), position3D.GetZ()));
    }

    else if (view == TPC_VIEW_W)
    {
        return CartesianVector(position3D.GetX(), 0.f, position3D.GetZ());
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::ProjectDirection(const Pandora &pandora, const CartesianVector &direction3D, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        return CartesianVector(direction3D.GetX(), 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->PYPZtoPU(direction3D.GetY(), direction3D.GetZ())).GetUnitVector();
    }

    else if (view == TPC_VIEW_V)
    {
        return CartesianVector(direction3D.GetX(), 0.f, LArGeometryHelper::GetLArTransformationPlugin(pandora)->PYPZtoPV(direction3D.GetY(), direction3D.GetZ())).GetUnitVector();
    }

    else if (view == TPC_VIEW_W)
    {
        return CartesianVector(direction3D.GetX(), 0.f, direction3D.GetZ()).GetUnitVector();
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::GetWireZPitch(const Pandora &pandora)
{
    return (LArGeometryHelper::GetLArPseudoLayerPlugin(pandora)->GetZPitch());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::GetWirePitch(const Pandora &pandora, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        return (LArGeometryHelper::GetLArPseudoLayerPlugin(pandora)->GetUPitch());
    }

    else if (view == TPC_VIEW_V)
    {
        return (LArGeometryHelper::GetLArPseudoLayerPlugin(pandora)->GetVPitch());
    }

    else if (view == TPC_VIEW_W)
    {
        return (LArGeometryHelper::GetLArPseudoLayerPlugin(pandora)->GetWPitch());
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::GetWireAxis(const Pandora &pandora, const HitType view)
{
    if (view == TPC_VIEW_U)
    {
        // CartesianVector(0.f, -m_sinU, m_cosU)
        return CartesianVector(0.f,
            LArGeometryHelper::GetLArTransformationPlugin(pandora)->PYPZtoPU(1.f, 0.f),
            LArGeometryHelper::GetLArTransformationPlugin(pandora)->PYPZtoPU(0.f, 1.f));
    }

    else if (view == TPC_VIEW_V)
    {
        // CartesianVector(0.f, +m_sinV, m_cosV)
        return CartesianVector(0.f,
        LArGeometryHelper::GetLArTransformationPlugin(pandora)->PYPZtoPV(1.f, 0.f),
        LArGeometryHelper::GetLArTransformationPlugin(pandora)->PYPZtoPV(0.f, 1.f));
    }

    else if (view == TPC_VIEW_W)
    {
        return CartesianVector(0.f, 0.f, 1.f);
    }

    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsInGap(const Pandora &pandora, const CartesianVector &testPoint2D, const HitType hitType, const float gapTolerance)
{
    // ATTN: input test point MUST be a 2D position vector
    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        // Check line gaps (need to input 2D hit type)
        const LineGap *const pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

        if (pLineGap)
        {
            if (pLineGap->IsInGap(testPoint2D, hitType, gapTolerance))
                return true;
        }

        // Check box gaps (need to input 3D hit type)
        const BoxGap *const pBoxGap = dynamic_cast<const BoxGap*>(pDetectorGap);

        if (pBoxGap)
        {
            if (pBoxGap->IsInGap(testPoint2D, TPC_3D, gapTolerance))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsInGap3D(const Pandora &pandora, const CartesianVector &testPoint3D, const HitType hitType, const float gapTolerance)
{
    const CartesianVector testPoint2D(LArGeometryHelper::ProjectPosition(pandora, testPoint3D, hitType));
    return LArGeometryHelper::IsInGap(pandora, testPoint2D, hitType, gapTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsXSamplingPointInGap(const Pandora &pandora, const float xSample, const TwoDSlidingFitResult &slidingFitResult,
    const float gapTolerance)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(slidingFitResult.GetCluster()));
    const CartesianVector minLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());

    const bool minLayerIsAtLowX(minLayerPosition.GetX() < maxLayerPosition.GetX());
    const CartesianVector &lowXCoordinate(minLayerIsAtLowX ? minLayerPosition : maxLayerPosition);
    const CartesianVector &highXCoordinate(minLayerIsAtLowX ? maxLayerPosition : minLayerPosition);

    if ((xSample > lowXCoordinate.GetX()) && (xSample < highXCoordinate.GetX()))
    {
        CartesianVector slidingFitPosition(0.f, 0.f, 0.f);

        if (STATUS_CODE_SUCCESS == slidingFitResult.GetGlobalFitPositionAtX(xSample, slidingFitPosition))
            return (LArGeometryHelper::IsInGap(pandora, slidingFitPosition, hitType, gapTolerance));
    }

    const CartesianVector lowXDirection(minLayerIsAtLowX ? slidingFitResult.GetGlobalMinLayerDirection() : slidingFitResult.GetGlobalMaxLayerDirection());
    const CartesianVector highXDirection(minLayerIsAtLowX ? slidingFitResult.GetGlobalMaxLayerDirection() : slidingFitResult.GetGlobalMinLayerDirection());

    const bool sampleIsNearerToLowX(std::fabs(xSample - lowXCoordinate.GetX()) < std::fabs(xSample - highXCoordinate.GetX()));
    const CartesianVector &startPosition(sampleIsNearerToLowX ? lowXCoordinate : highXCoordinate);
    const CartesianVector &startDirection(sampleIsNearerToLowX ? lowXDirection : highXDirection);

    if (std::fabs(startDirection.GetX()) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float pathLength((xSample - startPosition.GetX()) / startDirection.GetX());
    const CartesianVector samplingPoint(startPosition + startDirection * pathLength);

    return (LArGeometryHelper::IsInGap(pandora, samplingPoint, hitType, gapTolerance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArGeometryHelper::IsInGap(const Pandora &pandora, const CartesianVector &p1, const CartesianVector &p2, const HitType hitType,
    const float gapTolerance)
{
    // ATTN: input test points MUST be 2D position vectors
    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        // Check line gaps
        const LineGap *const pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

        if (pLineGap)
        {
            if (pLineGap->IsInGap(p1, hitType, gapTolerance) && pLineGap->IsInGap(p2, hitType, gapTolerance))
                return true;
        }

        // Check box gaps (need to input 3D hit type)
        const BoxGap *const pBoxGap = dynamic_cast<const BoxGap*>(pDetectorGap);

        if (pBoxGap)
        {
            if (pBoxGap->IsInGap(p1, TPC_3D, gapTolerance) && pBoxGap->IsInGap(p2, TPC_3D, gapTolerance))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateGapDeltaZ(const Pandora &pandora, const float minZ, const float maxZ, const HitType hitType)
{
    if (maxZ - minZ < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float gapFraction(0.f);

    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        // Check line gaps
        const LineGap *const pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

        if (pLineGap)
        {
            if (pLineGap->GetHitType() == hitType)
                gapFraction += LArGeometryHelper::CalculateOverlapFraction(pLineGap, minZ, maxZ);
        }
    }

    return gapFraction * (maxZ - minZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateGapDisplacement(const LineGap *const pLineGap, const float &z)
{
    const float deltaStartZ(z - pLineGap->GetLineStartZ());
    const float deltaEndZ(z - pLineGap->GetLineEndZ());

    return ((deltaStartZ < 0.f) ? std::fabs(deltaStartZ) : (deltaEndZ > 0.f) ? deltaEndZ : 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateGapDisplacement(const BoxGap *const pBoxGap, const CartesianVector &positionVector)
{
    bool inGap(false);
    float gapDisplacement(std::numeric_limits<float>::max());

    // This method should work for both 2D and 3D position vectors
    const CartesianVector relativePosition(positionVector - pBoxGap->GetVertex());

    // Calculate displacement from first pair of edges
    const float s1(pBoxGap->GetSide1().GetMagnitude());
    const float p1(relativePosition.GetDotProduct(pBoxGap->GetSide1().GetUnitVector()));
    const float x1((p1 < 0.f) ? std::fabs(p1) : (p1 > s1) ? (p1 - s1) : 0.f);

    if (x1 > std::numeric_limits<float>::epsilon())
    {
        inGap = true;
        gapDisplacement = std::min(gapDisplacement, x1);
    }

    // Calculate displacement from second pair of edges
    const float s2(pBoxGap->GetSide2().GetMagnitude());
    const float p2(relativePosition.GetDotProduct(pBoxGap->GetSide2().GetUnitVector()));
    const float x2((p2 < 0.f) ? std::fabs(p2) : (p2 > s2) ? (p2 - s2) : 0.f);

    if (x2 > std::numeric_limits<float>::epsilon())
    {
        inGap = true;
        gapDisplacement = std::min(gapDisplacement, x2);
    }

    // Calculate displacement from third pair of edges
    const float s3(pBoxGap->GetSide3().GetMagnitude());
    const float p3(relativePosition.GetDotProduct(pBoxGap->GetSide3().GetUnitVector()));
    const float x3((p3 < 0.f) ? std::fabs(p3) : (p3 > s3) ? (p3 - s3) : 0.f);

    if (x3 > std::numeric_limits<float>::epsilon())
    {
        inGap = true;
        gapDisplacement = std::min(gapDisplacement, x3);
    }

    if (inGap)
        return gapDisplacement;

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateOverlapFraction(const Pandora &pandora, const CartesianVector &p1, const CartesianVector &p2,
    const HitType hitType)
{
    float gapFraction(0.f);

    for (const DetectorGap *const pDetectorGap : pandora.GetGeometry()->GetDetectorGapList())
    {
        // Check line gaps
        const LineGap *const pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

        if (pLineGap)
        {
            if (pLineGap->GetHitType() == hitType)
                gapFraction += LArGeometryHelper::CalculateOverlapFraction(pLineGap, p1.GetZ(), p2.GetZ());
        }

        // Check box gaps
        const BoxGap *const pBoxGap = dynamic_cast<const BoxGap*>(pDetectorGap);

        if (pBoxGap)
        {
            gapFraction += LArGeometryHelper::CalculateOverlapFraction(pBoxGap, p1, p2);
        }
    }

    return std::min(gapFraction, 1.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateOverlapFraction(const LineGap *const pLineGap, const float &z1, const float &z2)
{
    const float pointMinZ(std::min(z1,z2));
    const float pointMaxZ(std::max(z1,z2));

    if (pLineGap->GetLineStartZ() > pointMaxZ || pLineGap->GetLineEndZ() < pointMinZ)
        return 0.f;

    const float overlapMinZ(std::max(pointMinZ,pLineGap->GetLineStartZ()));
    const float overlapMaxZ(std::min(pointMaxZ,pLineGap->GetLineEndZ()));

    return (overlapMaxZ - overlapMinZ) / (pointMaxZ - pointMinZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGeometryHelper::CalculateOverlapFraction(const BoxGap *const pBoxGap, const CartesianVector &p1, const CartesianVector &p2)
{
    // If input points are coincident, then just check whether they are inside the box gap
    if ((p2-p1).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
    {
        if (pBoxGap->IsInGap(p1, TPC_3D, 0.f) || pBoxGap->IsInGap(p2, TPC_3D, 0.f))
        {
            return 1.f;
        }
        else
        {
            return 0.f;
        }
    }

    // Cache properties of box gap to use in the calculations that follow
    const CartesianVector &s1(pBoxGap->GetSide1());
    const CartesianVector &s2(pBoxGap->GetSide2());
    const CartesianVector &s3(pBoxGap->GetSide3());
    const CartesianVector &v1(pBoxGap->GetVertex());

    const CartesianVector v2(v1 + s1 + s2 + s3);
    const CartesianVector traj((p2-p1).GetUnitVector());

    // Consider intersections with all six faces of the box gap
    float tmin(+std::numeric_limits<float>::max());
    float tmax(-std::numeric_limits<float>::max());

    try
    {
        const CartesianVector r1(LArGeometryHelper::CalculateIntersection(v1, s1, s2, p1, traj));
        const float t(traj.GetDotProduct((r1-p1)));
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    catch (StatusCodeException &) { }

    try
    {
        const CartesianVector r2(LArGeometryHelper::CalculateIntersection(v1, s2, s3, p1, traj));
        const float t(traj.GetDotProduct((r2-p1)));
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    catch (StatusCodeException &) { }

    try
    {
        const CartesianVector r3(LArGeometryHelper::CalculateIntersection(v1, s3, s1, p1, traj));
        const float t(traj.GetDotProduct((r3-p1)));
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    catch (StatusCodeException &) { }

    try
    {
        const CartesianVector r4(LArGeometryHelper::CalculateIntersection(v2, s1 * -1.f, s2 * -1.f, p1, traj));
        const float t(traj.GetDotProduct((r4-p1)));
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    catch (StatusCodeException &) { }

    try
    {
        const CartesianVector r5(LArGeometryHelper::CalculateIntersection(v2, s2 * -1.f, s3 * -1.f, p1, traj));
        const float t(traj.GetDotProduct((r5-p1)));
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    catch (StatusCodeException &) { }

    try
    {
        const CartesianVector r6(LArGeometryHelper::CalculateIntersection(v2, s3 * -1.f, s1 * -1.f, p1, traj));
        const float t(traj.GetDotProduct((r6-p1)));
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    catch (StatusCodeException &) { }

    // Calculate overlap between input pair of points and this box gap
    const float t0(0.f);
    const float t1((p2-p1).GetMagnitude());
    const float gmin(std::max(t0,tmin));
    const float gmax(std::min(t1,tmax));

    if (gmax - gmin > std::numeric_limits<float>::epsilon())
        return (gmax - gmin) / (t1 - t0);

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArGeometryHelper::CalculateIntersection(const CartesianVector &v0, const CartesianVector &s1, const CartesianVector &s2,
    const CartesianVector &r0, const CartesianVector &traj)
{
    // Calculate intersection of line with *finite* plane

    // Step (1): Calculate intersection of line with infinite plane
    //  Equation of line:    r = r0 + t * traj
    //  Equation of plane:  (r - v0).n = 0, where n = s1 x s2
    //  Intersection point:  t_int = (r0 - v0).n / traj.n
    //
    // Step (2): Check that line intersects the finite plane
    //  Not parallel:  traj.n != 0
    //  Finite plane:  0 < (r_int - v0).s1 < |s1|^2 and 0 < (r_int - v0).s2 < |s2|^2

    // Calculate intersection of line with the infinite plane
    const CartesianVector n(s1.GetCrossProduct(s2).GetUnitVector());
    const float n1(n.GetDotProduct(r0-v0));
    const float n2(n.GetDotProduct(traj));

    // Check that line and plane are not parallel
    if (std::fabs(n2) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Get the intersection point
    const CartesianVector p(r0 - traj * (n1/n2));

    // Check that line intersects with the finite plane
    const float t1((p - v0).GetDotProduct(s1.GetUnitVector()));
    const float t2((p - v0).GetDotProduct(s2.GetUnitVector()));

    if (t1 < 0.f || t2 < 0.f || t1 > s1.GetMagnitude() || t2 > s2.GetMagnitude())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return p;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArPseudoLayerPlugin *LArGeometryHelper::GetLArPseudoLayerPlugin(const Pandora &pandora)
{
    PseudoLayerInstanceMap::const_iterator iter = m_pseudolayerInstanceMap.find(&pandora);

    if (m_pseudolayerInstanceMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArTransformationPlugin *LArGeometryHelper::GetLArTransformationPlugin(const Pandora &pandora)
{
    TransformationInstanceMap::const_iterator iter = m_transformationInstanceMap.find(&pandora);

    if (m_transformationInstanceMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArGeometryHelper::SetLArPseudoLayerPlugin(const Pandora &pandora, const LArPseudoLayerPlugin *const pLArPseudoLayerPlugin)
{
    PseudoLayerInstanceMap::const_iterator iter = m_pseudolayerInstanceMap.find(&pandora);

    if (m_pseudolayerInstanceMap.end() != iter)
        return STATUS_CODE_ALREADY_INITIALIZED;

    if (!m_pseudolayerInstanceMap.insert(PseudoLayerInstanceMap::value_type(&pandora, pLArPseudoLayerPlugin)).second)
        return STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArGeometryHelper::SetLArTransformationPlugin(const Pandora &pandora, const LArTransformationPlugin *const pLArTransformationPlugin)
{
    TransformationInstanceMap::const_iterator iter = m_transformationInstanceMap.find(&pandora);

    if (m_transformationInstanceMap.end() != iter)
        return STATUS_CODE_ALREADY_INITIALIZED;

    if (!m_transformationInstanceMap.insert(TransformationInstanceMap::value_type(&pandora, pLArTransformationPlugin)).second)
        return STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
