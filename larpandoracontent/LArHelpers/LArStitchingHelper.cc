/**
 *  @file   larpandoracontent/LArHelpers/LArStitchingHelper.cc
 *
 *  @brief  Implementation of the helper class for multiple drift volumes
 *
 *  $Log: $
 */

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

const VolumeInfo &LArStitchingHelper::FindClosestVolume(const Pandora &pandora, const VolumeInfo &inputVolume, const bool checkPositive)
{
    const VolumeIdList &volumeIdList(MultiPandoraApi::GetVolumeIdList(&pandora));

    int closestId(-1);
    float closestSeparation(std::numeric_limits<float>::max());
    const float maxDisplacement(30.f); // TODO: 30cm should be fine, but can we do better than a hard-coded number here?

    for (const int volumeId : volumeIdList)
    {
        const VolumeInfo &checkVolume(MultiPandoraApi::GetVolumeInfo(&pandora, volumeId));

        if (inputVolume.GetIdNumber() == checkVolume.GetIdNumber())
            continue;

        if (checkPositive != (checkVolume.GetCenterX() > inputVolume.GetCenterX()))
            continue;

        const float deltaX(std::fabs(checkVolume.GetCenterX() - inputVolume.GetCenterX()));
        const float deltaY(std::fabs(checkVolume.GetCenterY() - inputVolume.GetCenterY()));
        const float deltaZ(std::fabs(checkVolume.GetCenterZ() - inputVolume.GetCenterZ()));

        if (deltaY > maxDisplacement || deltaZ > maxDisplacement)
            continue;

        if (deltaX < closestSeparation)
        {
            closestSeparation = deltaX;
            closestId = checkVolume.GetIdNumber();
        }
    }

    if (closestId < 0)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return MultiPandoraApi::GetVolumeInfo(&pandora, closestId);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::CanVolumesBeStitched(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume)
{
    // Require that drift volumes should have unique identifiers
    if (firstVolume.GetIdNumber() == secondVolume.GetIdNumber())
        return false;

    // ATTN: We assume that Pfos are crossing either an anode-anode boundary or a cathode-cathode boundary
    if (firstVolume.IsDriftInPositiveX() == secondVolume.IsDriftInPositiveX())
        return false;

    // Check if volumes are adjacent
    return LArStitchingHelper::AreVolumesAdjacent(firstVolume, secondVolume);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::AreVolumesAdjacent(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume)
{
    // Check the relative positions of the centres of each drift volume
    const float maxDisplacement(30.f); // TODO: 30cm should be fine, but can we do better than a hard-coded number here?
    const float widthX(0.5f * (firstVolume.GetWidthX() + secondVolume.GetWidthX()));
    const float deltaX(std::fabs(firstVolume.GetCenterX() - secondVolume.GetCenterX()));
    const float deltaY(std::fabs(firstVolume.GetCenterY() - secondVolume.GetCenterY()));
    const float deltaZ(std::fabs(firstVolume.GetCenterZ() - secondVolume.GetCenterZ()));

    if (std::fabs(deltaX-widthX) > maxDisplacement || deltaY > maxDisplacement || deltaZ > maxDisplacement)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArStitchingHelper::AreVolumesAdjacent(const Pandora &pandora, const VolumeInfo &firstVolume, const VolumeInfo &secondVolume)
{
    // Check if first volume is just upstream of second volume
    try
    {
        const VolumeInfo &firstVolumeCheck(LArStitchingHelper::FindClosestVolume(pandora, secondVolume, true));
        const VolumeInfo &secondVolumeCheck(LArStitchingHelper::FindClosestVolume(pandora, firstVolume, false));

        if (firstVolumeCheck.GetIdNumber() == firstVolume.GetIdNumber() && secondVolumeCheck.GetIdNumber() == secondVolume.GetIdNumber())
            return true;
    }
    catch (pandora::StatusCodeException& )
    {
    }

    // Check if second volume is just upstream of first volume
    try
    {
        const VolumeInfo &firstVolumeCheck(LArStitchingHelper::FindClosestVolume(pandora, secondVolume, false));
        const VolumeInfo &secondVolumeCheck(LArStitchingHelper::FindClosestVolume(pandora, firstVolume, true));

        if (firstVolumeCheck.GetIdNumber() == firstVolume.GetIdNumber() && secondVolumeCheck.GetIdNumber() == secondVolume.GetIdNumber())
            return true;
    }
    catch (pandora::StatusCodeException& )
    {
    }

    // Drift volumes aren't adjacent to each other
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetVolumeBoundaryCenterX(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume)
{
    if (!LArStitchingHelper::AreVolumesAdjacent(firstVolume, secondVolume))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (firstVolume.GetCenterX() < secondVolume.GetCenterX())
    {
        return 0.5 * ((firstVolume.GetCenterX() + 0.5 * firstVolume.GetWidthX()) +
                      (secondVolume.GetCenterX() - 0.5 * secondVolume.GetWidthX()));
    }
    else
    {
        return 0.5 * ((firstVolume.GetCenterX() - 0.5 * firstVolume.GetWidthX()) +
                      (secondVolume.GetCenterX() + 0.5 * secondVolume.GetWidthX()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetVolumeBoundaryWidthX(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume)
{
    if (!LArStitchingHelper::AreVolumesAdjacent(firstVolume, secondVolume))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (firstVolume.GetCenterX() < secondVolume.GetCenterX())
    {
        return ((secondVolume.GetCenterX() - 0.5 * secondVolume.GetWidthX()) -
                (firstVolume.GetCenterX() + 0.5 * firstVolume.GetWidthX()));
    }
    else
    {
        return ((firstVolume.GetCenterX() - 0.5 * firstVolume.GetWidthX()) -
                (secondVolume.GetCenterX() + 0.5 * secondVolume.GetWidthX()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::GetVolumeDisplacement(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume)
{
    const float deltaX(firstVolume.GetCenterX() - secondVolume.GetCenterX());
    const float deltaY(firstVolume.GetCenterY() - secondVolume.GetCenterY());
    const float deltaZ(firstVolume.GetCenterZ() - secondVolume.GetCenterZ());

    return std::sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArStitchingHelper::GetClosestVertices(const VolumeInfo &driftVolume1, const VolumeInfo &driftVolume2,
    const LArPointingCluster &pointingCluster1, const LArPointingCluster &pointingCluster2,
    LArPointingCluster::Vertex &closestVertex1, LArPointingCluster::Vertex &closestVertex2)
{
    // Check that drift volumes have different identifiers (just in case)
    if (driftVolume1.GetIdNumber() == driftVolume2.GetIdNumber())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Find the closest vertices based on relative X positions in drift volume
    const float dxVolume(driftVolume2.GetCenterX() - driftVolume1.GetCenterX());
    const float dx1(pointingCluster1.GetOuterVertex().GetPosition().GetX() - pointingCluster1.GetInnerVertex().GetPosition().GetX());
    const float dx2(pointingCluster2.GetOuterVertex().GetPosition().GetX() - pointingCluster2.GetInnerVertex().GetPosition().GetX());

    if (std::fabs(dx1) < std::numeric_limits<float>::epsilon() || std::fabs(dx2) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const bool useInner1((dxVolume > 0.f) == (dx1 < 0.f)); // xVol1 - OUTER - INNER - | - xVol2  [xVol2-xVol1>0; xOuter1-xInner1<0]
    const bool useInner2((dxVolume > 0.f) == (dx2 > 0.f)); // xVol1 - | - INNER - OUTER - xVol2  [xVol2-xVol1>0; xOuter2-xInner2>0]

    // Confirm that these really are the closest pair of vertices by checking the other possible pairs
    const LArPointingCluster::Vertex &nearVertex1(useInner1 ?  pointingCluster1.GetInnerVertex() : pointingCluster1.GetOuterVertex());
    const LArPointingCluster::Vertex &nearVertex2(useInner2 ?  pointingCluster2.GetInnerVertex() : pointingCluster2.GetOuterVertex());

    const LArPointingCluster::Vertex &farVertex1(useInner1 ?  pointingCluster1.GetOuterVertex() : pointingCluster1.GetInnerVertex());
    const LArPointingCluster::Vertex &farVertex2(useInner2 ?  pointingCluster2.GetOuterVertex() : pointingCluster2.GetInnerVertex());

    const float dxNearNear(0.f);
    const float dyNearNear(nearVertex1.GetPosition().GetY() - nearVertex2.GetPosition().GetY());
    const float dzNearNear(nearVertex1.GetPosition().GetZ() - nearVertex2.GetPosition().GetZ());
    const float drNearNearSquared(dxNearNear * dxNearNear + dyNearNear * dyNearNear + dzNearNear * dzNearNear);

    const float dxNearFar(std::fabs(dx2));
    const float dyNearFar(nearVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());
    const float dzNearFar(nearVertex1.GetPosition().GetZ() - farVertex2.GetPosition().GetZ());
    const float drNearFarSquared(dxNearFar * dxNearFar + dyNearFar * dyNearFar + dzNearFar * dzNearFar);

    const float dxFarNear(std::fabs(dx1));
    const float dyFarNear(farVertex1.GetPosition().GetY() - nearVertex2.GetPosition().GetY());
    const float dzFarNear(farVertex1.GetPosition().GetZ() - nearVertex2.GetPosition().GetZ());
    const float drFarNearSquared(dxFarNear * dxFarNear + dyFarNear * dyFarNear + dzFarNear * dzFarNear);

    const float dxFarFar(std::fabs(dx1) + std::fabs(dx2));
    const float dyFarFar(farVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());
    const float dzFarFar(farVertex1.GetPosition().GetZ() - farVertex2.GetPosition().GetZ());
    const float drFarFarSquared(dxFarFar * dxFarFar + dyFarFar * dyFarFar + dzFarFar * dzFarFar);

    if (drNearNearSquared > std::min(drFarFarSquared, std::min(drNearFarSquared, drFarNearSquared)))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    closestVertex1 = nearVertex1;
    closestVertex2 = nearVertex2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArStitchingHelper::CalculateX0(const VolumeInfo &firstVolume, const VolumeInfo &secondVolume,
    const LArPointingCluster::Vertex &firstVertex, const LArPointingCluster::Vertex &secondVertex)
{
    // Require that drift volumes should have unique identifiers
    if (firstVolume.GetIdNumber() == secondVolume.GetIdNumber())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // ATTN: Assume that Pfos are crossing either an anode-anode boundary or a cathode-cathode boundary
    if (firstVolume.IsDriftInPositiveX() == secondVolume.IsDriftInPositiveX())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Assume that Pfos have opposite direction components in x, and have some direction component in the y-z plane
    const CartesianVector firstDirectionX(firstVertex.GetDirection().GetX(), 0.f, 0.f);
    const CartesianVector secondDirectionX(secondVertex.GetDirection().GetX(), 0.f, 0.f);

    if (std::fabs(firstDirectionX.GetX()) > 1.f - std::numeric_limits<float>::epsilon() ||
        std::fabs(secondDirectionX.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (-firstDirectionX.GetDotProduct(secondDirectionX) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Calculate displacement in x from relative displacement in y-z plane
    const CartesianVector firstPositionYZ(0.f, firstVertex.GetPosition().GetY(), firstVertex.GetPosition().GetZ());
    const CartesianVector firstDirectionYZ(0.f, firstVertex.GetDirection().GetY(), firstVertex.GetDirection().GetZ());

    const CartesianVector secondPositionYZ(0.f, secondVertex.GetPosition().GetY(), secondVertex.GetPosition().GetZ());
    const CartesianVector secondDirectionYZ(0.f, secondVertex.GetDirection().GetY(), secondVertex.GetDirection().GetZ());

    const float firstDirectionYZmag(firstDirectionYZ.GetMagnitude());
    const float secondDirectionYZmag(secondDirectionYZ.GetMagnitude());

    if (firstDirectionYZmag < std::numeric_limits<float>::epsilon() || secondDirectionYZmag < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float R1(firstDirectionYZ.GetUnitVector().GetDotProduct(firstPositionYZ - secondPositionYZ) / firstDirectionYZmag);
    const float X1(-1.f * firstDirectionX.GetUnitVector().GetDotProduct(secondVertex.GetPosition() - (firstVertex.GetPosition() - firstVertex.GetDirection() * R1)));

    const float R2(secondDirectionYZ.GetUnitVector().GetDotProduct(secondPositionYZ - firstPositionYZ) / secondDirectionYZmag);
    const float X2(-1.f * secondDirectionX.GetUnitVector().GetDotProduct(firstVertex.GetPosition() - (secondVertex.GetPosition() - secondVertex.GetDirection() * R2)));

    // ATTN: By convention, X0 is half the displacement in x (because both Pfos will be corrected)
    return (X1 + X2) * 0.25f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArStitchingHelper::GetCorrectedPosition(const VolumeInfo &driftVolume, const float x0, const CartesianVector &inputPosition)
{
  return (inputPosition + CartesianVector(driftVolume.IsDriftInPositiveX() ? -x0 : x0, 0.f, 0.f));
}

} // namespace lar_content
