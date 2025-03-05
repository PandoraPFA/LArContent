/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLLaterTierHierarchyTool.cc
 *
 *  @brief  Implementation of the DL later tier hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLLaterTierHierarchyTool.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::DLLaterTierNetworkParams::AddCommonParamsToInput(int &insertIndex, LArDLHelper::TorchInput &modelInput) const
{
    for (const float param : {m_parentTrackScore, m_childTrackScore, m_parentNSpacepoints, m_childNSpacepoints, m_separation3D})
    {
        modelInput[0][insertIndex] = param;
        ++insertIndex;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::DLLaterTierNetworkParams::AddOrientationParamsToInput(int &insertIndex, LArDLHelper::TorchInput &modelInput) const
{
    for (const float param : {m_parentNuVertexSep, m_childNuVertexSep, m_parentEndRegionNHits, m_parentEndRegionNParticles,
             m_parentEndRegionRToWall, m_vertexSeparation, m_doesChildConnect, m_overshootStartDCA, m_overshootStartL, m_overshootEndDCA,
             m_overshootEndL, m_childCPDCA, m_childCPExtrapDistance, m_childCPLRatio, m_parentCPNUpstreamHits, m_parentCPNDownstreamHits,
             m_parentCPNHitRatio, m_parentCPEigenvalueRatio, m_parentCPOpeningAngle, m_parentIsPOIClosestToNu, m_childIsPOIClosestToNu})
    {
        modelInput[0][insertIndex] = param;
        ++insertIndex;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DLLaterTierHierarchyTool::DLLaterTierHierarchyTool() :
    DLBaseHierarchyTool(),
    m_trainingMode(false),
    m_trajectoryStepSize(1.f),
    m_connectionBuffer(50.f),
    m_searchRegion(10.f),
    m_normalise(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLLaterTierHierarchyTool::Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pNeutrinoPfo,
    const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo,
    std::vector<DLLaterTierNetworkParams> &networkParamVector, float &laterTierScore)
{
    networkParamVector.clear();
    laterTierScore = std::numeric_limits<float>::lowest();

    this->SetDetectorBoundaries();

    // Get the common params i.e. independent of postulated orientation
    const std::pair<float, float> trackScoreParams(
        {LArPfoHelper::GetTrackScore(parentHierarchyPfo.GetPfo()), LArPfoHelper::GetTrackScore(childHierarchyPfo.GetPfo())});
    const std::pair<float, float> nSpacepointsParams(
        {LArPfoHelper::GetNumberOfThreeDHits(parentHierarchyPfo.GetPfo()), LArPfoHelper::GetNumberOfThreeDHits(childHierarchyPfo.GetPfo())});
    const float separation3D(LArPfoHelper::GetThreeDSeparation(parentHierarchyPfo.GetPfo(), childHierarchyPfo.GetPfo()));

    const bool isChildTrack(childHierarchyPfo.GetPfo()->GetParticleId() == 13);
    std::vector<bool> childOrientationVector(
        isChildTrack ? std::vector<bool>({true, false}) : std::vector<bool>({this->IsShowerVertexUpstream(parentHierarchyPfo, childHierarchyPfo)}));

    // Set network params
    for (const bool &useUpstreamForParent : {true, false})
    {
        for (const bool &useUpstreamForChild : childOrientationVector)
        {
            DLLaterTierNetworkParams edgeParams;

            // Set common params
            this->SetCommonParams(trackScoreParams, nSpacepointsParams, separation3D, edgeParams);

            // Get orientation dependent params
            const StatusCode statusCode(this->CalculateNetworkVariables(
                pAlgorithm, parentHierarchyPfo, childHierarchyPfo, pNeutrinoPfo, useUpstreamForParent, useUpstreamForChild, edgeParams));

            if (statusCode != STATUS_CODE_SUCCESS)
                return statusCode;

            // Add these to our output vector
            networkParamVector.emplace_back(edgeParams);
        }
    }

    // Now run the model!
    if (!m_trainingMode)
    {
        if (isChildTrack)
            laterTierScore =
                this->ClassifyTrackTrack(networkParamVector.at(0), networkParamVector.at(1), networkParamVector.at(2), networkParamVector.at(3));
        else
            laterTierScore = this->ClassifyTrackShower(networkParamVector.at(0), networkParamVector.at(1));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLLaterTierHierarchyTool::IsShowerVertexUpstream(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const
{
    const bool isUpstreamSet(childHierarchyPfo.GetUpstreamPoint().IsSet());
    const bool isDownstreamSet(childHierarchyPfo.GetDownstreamPoint().IsSet());

    if (isUpstreamSet && !isDownstreamSet)
        return true;

    if (!isUpstreamSet && isDownstreamSet)
        return false;

    // By design, showers should have at least one set
    if (!isUpstreamSet && !isDownstreamSet)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    CartesianPointVector parentPositions3D;
    LArPfoHelper::GetCoordinateVector(parentHierarchyPfo.GetPfo(), TPC_3D, parentPositions3D);

    float upstreamSepSq(std::numeric_limits<float>::max());
    float downstreamSepSq(std::numeric_limits<float>::max());

    for (const CartesianVector &parentPosition3D : parentPositions3D)
    {
        const float thisUpstreamSepSq((parentPosition3D - childHierarchyPfo.GetUpstreamPoint().GetPosition()).GetMagnitudeSquared());
        const float thisDownstreamSepSq((parentPosition3D - childHierarchyPfo.GetDownstreamPoint().GetPosition()).GetMagnitudeSquared());

        upstreamSepSq = std::min(upstreamSepSq, thisUpstreamSepSq);
        downstreamSepSq = std::min(downstreamSepSq, thisDownstreamSepSq);
    }

    return (upstreamSepSq < downstreamSepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetCommonParams(const std::pair<float, float> &trackScoreParams,
    const std::pair<float, float> &nSpacepointsParams, const float separation3D, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    laterTierNetworkParams.m_parentTrackScore = trackScoreParams.first;
    laterTierNetworkParams.m_childTrackScore = trackScoreParams.second;
    laterTierNetworkParams.m_parentNSpacepoints = nSpacepointsParams.first;
    laterTierNetworkParams.m_childNSpacepoints = nSpacepointsParams.second;
    laterTierNetworkParams.m_separation3D = separation3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLLaterTierHierarchyTool::CalculateNetworkVariables(const Algorithm *const pAlgorithm, const HierarchyPfo &parentHierarchyPfo,
    const HierarchyPfo &childHierarchyPfo, const ParticleFlowObject *const pNeutrinoPfo, const bool useUpstreamForParent,
    const bool useUpstreamForChild, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    // Pick out neutrino vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;

    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(
        pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    // Pick out the correct pfo vertex/direction
    // Parent POI is a postulate endpoint
    // Child POI is a postulate startpoint
    const ExtremalPoint &parentStart(useUpstreamForParent ? parentHierarchyPfo.GetDownstreamPoint() : parentHierarchyPfo.GetUpstreamPoint());
    const ExtremalPoint &parentEnd(useUpstreamForParent ? parentHierarchyPfo.GetUpstreamPoint() : parentHierarchyPfo.GetDownstreamPoint());
    const ExtremalPoint &childStart(useUpstreamForChild ? childHierarchyPfo.GetUpstreamPoint() : childHierarchyPfo.GetDownstreamPoint());

    // Check that we could actually calculate pfo vertex/direction, return if not
    if (!parentStart.IsSet())
    {
        return STATUS_CODE_NOT_INITIALIZED;
    }
    if (!parentEnd.IsSet())
    {
        return STATUS_CODE_NOT_INITIALIZED;
    }
    if (!childStart.IsSet())
    {
        return STATUS_CODE_NOT_INITIALIZED;
    }

    // Set laterTierNetworkParams
    laterTierNetworkParams.m_parentIsPOIClosestToNu = useUpstreamForParent ? 1.f : 0.f;
    laterTierNetworkParams.m_childIsPOIClosestToNu = useUpstreamForChild ? 1.f : 0.f;
    this->SetVertexParams(nuVertex, parentStart.GetPosition(), parentEnd.GetPosition(), childStart.GetPosition(), laterTierNetworkParams);
    this->SetEndRegionParams(pAlgorithm, parentHierarchyPfo.GetPfo(), parentEnd.GetPosition(), laterTierNetworkParams);
    this->SetConnectionParams(parentHierarchyPfo, parentStart.GetPosition(), childStart, laterTierNetworkParams);
    this->SetOvershootParams(parentStart, parentEnd, childStart, laterTierNetworkParams);
    this->SetParentConnectionPointVars(parentHierarchyPfo, laterTierNetworkParams);

    // Normalise
    if (m_normalise)
        this->NormaliseNetworkParams(laterTierNetworkParams);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetVertexParams(const CartesianVector &nuVertex, const CartesianVector &parentStartPos,
    const CartesianVector &parentEndPos, const CartesianVector &childStartPos, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    laterTierNetworkParams.m_parentNuVertexSep = (parentStartPos - nuVertex).GetMagnitude();
    laterTierNetworkParams.m_childNuVertexSep = (childStartPos - nuVertex).GetMagnitude();
    laterTierNetworkParams.m_vertexSeparation = (parentEndPos - childStartPos).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetEndRegionParams(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pParentPfo,
    const CartesianVector &parentEndPos, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    std::pair<float, float> endRegionParams(this->GetParticleInfoAboutPfoPosition(pAlgorithm, pParentPfo, parentEndPos));

    laterTierNetworkParams.m_parentEndRegionNHits = endRegionParams.first;
    laterTierNetworkParams.m_parentEndRegionNParticles = endRegionParams.second;

    this->SetEndRegionRToWall(parentEndPos, laterTierNetworkParams);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetEndRegionRToWall(const CartesianVector &parentEndPos, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    float closestDistance(std::numeric_limits<float>::max());

    for (int iAxis : {0, 1, 2})
    {
        const float parentCoord(iAxis == 0 ? parentEndPos.GetX() : iAxis == 1 ? parentEndPos.GetY() : parentEndPos.GetZ());
        const std::pair<float, float> &boundary(iAxis == 0 ? m_detectorBoundaries.m_xBoundaries
                : iAxis == 1                               ? m_detectorBoundaries.m_yBoundaries
                                                           : m_detectorBoundaries.m_zBoundaries);
        const float dimensionMin(std::min((parentCoord - boundary.first), (boundary.second - parentCoord)));

        closestDistance = std::min(closestDistance, dimensionMin);
    }

    laterTierNetworkParams.m_parentEndRegionRToWall = closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetConnectionParams(const HierarchyPfo &parentHierarchyPfo, const CartesianVector &parentStartPos,
    const ExtremalPoint &childStart, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    // Get fit of parent
    const ThreeDSlidingFitResult &slidingFitResult(parentHierarchyPfo.GetSlidingFitResult());

    // Now walk along track, starting from parentStartPos
    // I want this to be interpreted as l = 0 so have to alter pandora def
    const bool startFromMin((parentStartPos - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitudeSquared() <
        (parentStartPos - slidingFitResult.GetGlobalMaxLayerPosition()).GetMagnitudeSquared());
    const float minL(slidingFitResult.GetLongitudinalDisplacement(slidingFitResult.GetGlobalMinLayerPosition()));
    const float maxL(slidingFitResult.GetLongitudinalDisplacement(slidingFitResult.GetGlobalMaxLayerPosition()));

    if (std::fabs(minL - maxL) < std::numeric_limits<float>::epsilon())
        return;

    const int nSteps(std::floor(std::fabs(maxL - minL) / m_trajectoryStepSize));

    // Keep track of these in loop
    float slidingL(startFromMin ? minL : 0.f);
    float trackLength(0.f);
    float minDCA(std::numeric_limits<float>::max());
    float lengthAtConnection(0.f);

    laterTierNetworkParams.m_doesChildConnect = 0.f;

    for (int i = 0; i < nSteps; i++)
    {
        float firstL(slidingL), secondL(slidingL + m_trajectoryStepSize);
        float firstEvalL(startFromMin ? firstL : (maxL - firstL));
        float secondEvalL(startFromMin ? secondL : (maxL - secondL));
        CartesianVector firstPosition(0.f, 0.f, 0.f);
        CartesianVector secondPosition(0.f, 0.f, 0.f);

        if ((slidingFitResult.GetGlobalFitPosition(firstEvalL, firstPosition) != STATUS_CODE_SUCCESS) ||
            (slidingFitResult.GetGlobalFitPosition(secondEvalL, secondPosition) != STATUS_CODE_SUCCESS))
        {
            slidingL += m_trajectoryStepSize;
            continue;
        }

        // Get trajectory position
        const CartesianVector midPoint((secondPosition + firstPosition) * 0.5f);

        // Extrap child vertex to parent position
        const std::pair<CartesianVector, bool> extrap(this->ExtrapolateChildToParent(midPoint, childStart));
        const float thisDCA((midPoint - extrap.first).GetMagnitude());

        if (thisDCA < minDCA)
        {
            if (this->DoesConnect(firstPosition, secondPosition, extrap.first))
            {
                minDCA = thisDCA;
                laterTierNetworkParams.m_doesChildConnect = 1.f;
                laterTierNetworkParams.m_connectionPoint = midPoint;
                laterTierNetworkParams.m_connectionDirection = (secondPosition - firstPosition).GetUnitVector();
                laterTierNetworkParams.m_childCPDCA = thisDCA;
                laterTierNetworkParams.m_childCPExtrapDistance =
                    (extrap.first - childStart.GetPosition()).GetMagnitude() * (extrap.second ? 1.f : (-1.f));
                lengthAtConnection = trackLength + (midPoint - firstPosition).GetMagnitude();
            }
        }

        slidingL += m_trajectoryStepSize;
        trackLength += (secondPosition - firstPosition).GetMagnitude();
    }

    if (std::fabs(laterTierNetworkParams.m_doesChildConnect - 1.f) < std::numeric_limits<float>::epsilon())
        laterTierNetworkParams.m_childCPLRatio = lengthAtConnection / trackLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<CartesianVector, bool> DLLaterTierHierarchyTool::ExtrapolateChildToParent(const CartesianVector &parentPos, const ExtremalPoint &childStart) const
{
    const double extrapDistance((parentPos - childStart.GetPosition()).GetDotProduct(childStart.GetDirection()));
    const bool backwardsExtrap(extrapDistance < 0.f);
    const CartesianVector extrapPoint(childStart.GetPosition() + (childStart.GetDirection() * extrapDistance));

    return std::pair<CartesianVector, bool>({extrapPoint, backwardsExtrap});
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLLaterTierHierarchyTool::DoesConnect(const CartesianVector &boundary1, const CartesianVector &boundary2, const CartesianVector &testPoint) const
{
    const CartesianVector displacement(boundary2 - boundary1);
    const double segmentLength(displacement.GetMagnitude());
    const CartesianVector unitVector(displacement * (1.f / segmentLength));
    const double l(unitVector.GetDotProduct(testPoint - boundary1));

    if ((l < 0.f) || (l > segmentLength))
        return false;

    const double t(unitVector.GetCrossProduct(testPoint - boundary1).GetMagnitudeSquared());

    if (t > (m_connectionBuffer * m_connectionBuffer))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetOvershootParams(const ExtremalPoint &parentStart, const ExtremalPoint &parentEnd,
    const ExtremalPoint &childStart, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    if (std::fabs(laterTierNetworkParams.m_doesChildConnect - 1.f) < std::numeric_limits<float>::epsilon())
        return;

    // Let's start at the very beginning, a very good place to start
    const std::pair<CartesianVector, bool> extrapStart(this->ExtrapolateChildToParent(parentStart.GetPosition(), childStart));
    laterTierNetworkParams.m_overshootStartDCA =
        (extrapStart.first - parentStart.GetPosition()).GetMagnitude() * (extrapStart.second ? 1.f : (-1.f));
    laterTierNetworkParams.m_overshootStartL =
        std::fabs((childStart.GetPosition() - parentStart.GetPosition()).GetDotProduct(parentStart.GetDirection()));

    // Now end
    const std::pair<CartesianVector, bool> extrapEnd(this->ExtrapolateChildToParent(parentEnd.GetPosition(), childStart));
    laterTierNetworkParams.m_overshootEndDCA = (extrapEnd.first - parentEnd.GetPosition()).GetMagnitude() * (extrapEnd.second ? 1.f : (-1.f));
    laterTierNetworkParams.m_overshootEndL =
        std::fabs((childStart.GetPosition() - parentEnd.GetPosition()).GetDotProduct(parentEnd.GetDirection()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::SetParentConnectionPointVars(const HierarchyPfo &parentHierarchyPfo, DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    if (std::fabs(laterTierNetworkParams.m_doesChildConnect - 1.f) > std::numeric_limits<float>::epsilon())
        return;

    // Use direction to current direction to split spacepoints into two groups
    CartesianPointVector parentPositions3D;
    LArPfoHelper::GetCoordinateVector(parentHierarchyPfo.GetPfo(), TPC_3D, parentPositions3D);

    CartesianPointVector upstreamGroup, downstreamGroup;

    for (const CartesianVector &parentPosition : parentPositions3D)
    {
        const CartesianVector displacement(parentPosition - laterTierNetworkParams.m_connectionPoint);

        // Do box as quicker than circle with mag.
        if ((std::fabs(displacement.GetX()) > m_searchRegion) || (std::fabs(displacement.GetY()) > m_searchRegion) ||
            (std::fabs(displacement.GetZ()) > m_searchRegion))
        {
            continue;
        }

        // Assign to group
        const float thisL(laterTierNetworkParams.m_connectionDirection.GetDotProduct(displacement));
        CartesianPointVector &group(thisL > 0.f ? downstreamGroup : upstreamGroup);
        group.emplace_back(parentPosition);
    }

    laterTierNetworkParams.m_parentCPNUpstreamHits = upstreamGroup.size();
    laterTierNetworkParams.m_parentCPNDownstreamHits = downstreamGroup.size();

    if (upstreamGroup.empty() || downstreamGroup.empty())
        return;

    laterTierNetworkParams.m_parentCPNHitRatio = laterTierNetworkParams.m_parentCPNDownstreamHits / laterTierNetworkParams.m_parentCPNUpstreamHits;

    // Now PCA magic
    try
    {
        CartesianVector centroidUp(0.f, 0.f, 0.f);
        CartesianVector centroidDown(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecsUp, eigenVecsDown;
        LArPcaHelper::EigenValues eigenValuesUp(0.f, 0.f, 0.f);
        LArPcaHelper::EigenValues eigenValuesDown(0.f, 0.f, 0.f);

        LArPcaHelper::RunPca(upstreamGroup, centroidUp, eigenValuesUp, eigenVecsUp);
        LArPcaHelper::RunPca(downstreamGroup, centroidDown, eigenValuesDown, eigenVecsDown);

        // Get opening angle from first eigenvectors (this is the longitudinal one) - straight would be around 180
        laterTierNetworkParams.m_parentCPOpeningAngle = eigenVecsUp.at(0).GetOpeningAngle(eigenVecsDown.at(0)) * (180.f / 3.14);

        // Get average transverse eigenvalues, get ratio
        const float avTransEVaUp((eigenValuesUp.GetY() + eigenValuesUp.GetZ()) * 0.5f);
        const float avTransEVaDown((eigenValuesDown.GetY() + eigenValuesDown.GetZ()) * 0.5f);

        laterTierNetworkParams.m_parentCPEigenvalueRatio = avTransEVaUp < std::numeric_limits<double>::epsilon() ? 0.f : avTransEVaDown / avTransEVaUp;
    }
    catch (...)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLLaterTierHierarchyTool::NormaliseNetworkParams(DLLaterTierNetworkParams &laterTierNetworkParams) const
{
    this->NormaliseNetworkParam(m_normLimits.m_trackScoreMin, m_normLimits.m_trackScoreMax, laterTierNetworkParams.m_parentTrackScore);
    this->NormaliseNetworkParam(m_normLimits.m_trackScoreMin, m_normLimits.m_trackScoreMax, laterTierNetworkParams.m_childTrackScore);
    this->NormaliseNetworkParam(m_normLimits.m_nSpacepointsMin, m_normLimits.m_nSpacepointsMax, laterTierNetworkParams.m_parentNSpacepoints);
    this->NormaliseNetworkParam(m_normLimits.m_nSpacepointsMin, m_normLimits.m_nSpacepointsMax, laterTierNetworkParams.m_childNSpacepoints);
    this->NormaliseNetworkParam(m_normLimits.m_separation3DMin, m_normLimits.m_separation3DMax, laterTierNetworkParams.m_separation3D);
    this->NormaliseNetworkParam(m_normLimits.m_nuVertexSepMin, m_normLimits.m_nuVertexSepMax, laterTierNetworkParams.m_parentNuVertexSep);
    this->NormaliseNetworkParam(m_normLimits.m_nuVertexSepMin, m_normLimits.m_nuVertexSepMax, laterTierNetworkParams.m_childNuVertexSep);
    this->NormaliseNetworkParam(m_normLimits.m_parentEndRegionNHitsMin, m_normLimits.m_parentEndRegionNHitsMax, laterTierNetworkParams.m_parentEndRegionNHits);
    this->NormaliseNetworkParam(m_normLimits.m_parentEndRegionNParticlesMin, m_normLimits.m_parentEndRegionNParticlesMax,
        laterTierNetworkParams.m_parentEndRegionNParticles);
    this->NormaliseNetworkParam(
        m_normLimits.m_parentEndRegionRToWallMin, m_normLimits.m_parentEndRegionRToWallMax, laterTierNetworkParams.m_parentEndRegionRToWall);
    this->NormaliseNetworkParam(m_normLimits.m_vertexSepMin, m_normLimits.m_vertexSepMax, laterTierNetworkParams.m_vertexSeparation);
    this->NormaliseNetworkParam(m_normLimits.m_doesChildConnectMin, m_normLimits.m_doesChildConnectMax, laterTierNetworkParams.m_doesChildConnect);
    this->NormaliseNetworkParam(m_normLimits.m_overshootDCAMin, m_normLimits.m_overshootDCAMax, laterTierNetworkParams.m_overshootStartDCA);
    this->NormaliseNetworkParam(m_normLimits.m_overshootLMin, m_normLimits.m_overshootLMax, laterTierNetworkParams.m_overshootStartL);
    this->NormaliseNetworkParam(m_normLimits.m_overshootDCAMin, m_normLimits.m_overshootDCAMax, laterTierNetworkParams.m_overshootEndDCA);
    this->NormaliseNetworkParam(m_normLimits.m_overshootLMin, m_normLimits.m_overshootLMax, laterTierNetworkParams.m_overshootEndL);
    this->NormaliseNetworkParam(m_normLimits.m_childCPDCAMin, m_normLimits.m_childCPDCAMax, laterTierNetworkParams.m_childCPDCA);
    this->NormaliseNetworkParam(
        m_normLimits.m_childCPExtrapDistanceMin, m_normLimits.m_childCPExtrapDistanceMax, laterTierNetworkParams.m_childCPExtrapDistance);
    this->NormaliseNetworkParam(m_normLimits.m_childCPLRatioMin, m_normLimits.m_childCPLRatioMax, laterTierNetworkParams.m_childCPLRatio);
    this->NormaliseNetworkParam(m_normLimits.m_parentCPNHitsMin, m_normLimits.m_parentCPNHitsMax, laterTierNetworkParams.m_parentCPNUpstreamHits);
    this->NormaliseNetworkParam(m_normLimits.m_parentCPNHitsMin, m_normLimits.m_parentCPNHitsMax, laterTierNetworkParams.m_parentCPNDownstreamHits);
    this->NormaliseNetworkParam(m_normLimits.m_parentCPNHitRatioMin, m_normLimits.m_parentCPNHitRatioMax, laterTierNetworkParams.m_parentCPNHitRatio);
    this->NormaliseNetworkParam(m_normLimits.m_parentCPEigenvalueRatioMin, m_normLimits.m_parentCPEigenvalueRatioMax,
        laterTierNetworkParams.m_parentCPEigenvalueRatio);
    this->NormaliseNetworkParam(m_normLimits.m_parentCPOpeningAngleMin, m_normLimits.m_parentCPOpeningAngleMax, laterTierNetworkParams.m_parentCPOpeningAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLLaterTierHierarchyTool::ClassifyTrackTrack(const DLLaterTierNetworkParams &edgeParamsUpUp, const DLLaterTierNetworkParams &edgeParamsUpDown,
    const DLLaterTierNetworkParams &edgeParamsDownUp, const DLLaterTierNetworkParams &edgeParamsDownDown)
{
    // Invoke branch model for each edge
    const FloatVector outputUpUp(this->ClassifyTrackTrackEdge(edgeParamsUpUp, edgeParamsUpDown, edgeParamsDownUp, edgeParamsDownDown));
    const FloatVector outputUpDown(this->ClassifyTrackTrackEdge(edgeParamsUpDown, edgeParamsDownUp, edgeParamsDownDown, edgeParamsUpUp));
    const FloatVector outputDownUp(this->ClassifyTrackTrackEdge(edgeParamsDownUp, edgeParamsDownDown, edgeParamsUpUp, edgeParamsUpDown));
    const FloatVector outputDownDown(this->ClassifyTrackTrackEdge(edgeParamsDownDown, edgeParamsUpUp, edgeParamsUpDown, edgeParamsDownUp));

    // Invoke classifier model for final output
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 12}, input);

    int insertIndex(0);

    for (const FloatVector &edgeOutput : {outputUpUp, outputUpDown, outputDownUp, outputDownDown})
    {
        for (int i = 0; i < 3; ++i)
        {
            input[0][insertIndex] = edgeOutput.at(i);
            ++insertIndex;
        }
    }

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_trackTrackClassifierModel, {input}, output);
    const torch::TensorAccessor<float, 2> outputAccessor = output.accessor<float, 2>();

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector DLLaterTierHierarchyTool::ClassifyTrackTrackEdge(const DLLaterTierNetworkParams &edgeParams,
    const DLLaterTierNetworkParams &otherEdgeParams1, const DLLaterTierNetworkParams &otherEdgeParams2, const DLLaterTierNetworkParams &otherEdgeParams3)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 89}, input);

    int insertIndex(0);
    edgeParams.AddCommonParamsToInput(insertIndex, input);
    edgeParams.AddOrientationParamsToInput(insertIndex, input);
    otherEdgeParams1.AddOrientationParamsToInput(insertIndex, input);
    otherEdgeParams2.AddOrientationParamsToInput(insertIndex, input);
    otherEdgeParams3.AddOrientationParamsToInput(insertIndex, input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_trackTrackBranchModel, {input}, output);

    const torch::TensorAccessor<float, 2> outputAccessor(output.accessor<float, 2>());

    return {outputAccessor[0][0], outputAccessor[0][1], outputAccessor[0][2]};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLLaterTierHierarchyTool::ClassifyTrackShower(const DLLaterTierNetworkParams &edgeParamsUp, const DLLaterTierNetworkParams &edgeParamsDown)
{
    // Invoke branch model for each edge
    const FloatVector outputUp(this->ClassifyTrackShowerEdge(edgeParamsUp, edgeParamsDown));
    const FloatVector outputDown(this->ClassifyTrackShowerEdge(edgeParamsDown, edgeParamsUp));

    // Invoke classifier model for final output
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 6}, input);

    int insertIndex(0);

    for (const FloatVector &edgeOutput : {outputUp, outputDown})
    {
        for (int i = 0; i < 3; ++i)
        {
            input[0][insertIndex] = edgeOutput.at(i);
            ++insertIndex;
        }
    }

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_trackShowerClassifierModel, {input}, output);
    const torch::TensorAccessor<float, 2> outputAccessor = output.accessor<float, 2>();

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector DLLaterTierHierarchyTool::ClassifyTrackShowerEdge(const DLLaterTierNetworkParams &edgeParams, const DLLaterTierNetworkParams &otherEdgeParams)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 47}, input);

    int insertIndex(0);
    edgeParams.AddCommonParamsToInput(insertIndex, input);
    edgeParams.AddOrientationParamsToInput(insertIndex, input);
    otherEdgeParams.AddOrientationParamsToInput(insertIndex, input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_trackShowerBranchModel, {input}, output);

    const torch::TensorAccessor<float, 2> outputAccessor(output.accessor<float, 2>());

    return {outputAccessor[0][0], outputAccessor[0][1], outputAccessor[0][2]};
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLLaterTierHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrajectoryStepSize", m_trajectoryStepSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConnectionBuffer", m_connectionBuffer));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SearchRegion", m_searchRegion));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Normalise", m_normalise));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackScoreMin", m_normLimits.m_trackScoreMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackScoreMax", m_normLimits.m_trackScoreMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMin", m_normLimits.m_nSpacepointsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMax", m_normLimits.m_nSpacepointsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Separation3DMin", m_normLimits.m_separation3DMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Separation3DMax", m_normLimits.m_separation3DMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexSepMin", m_normLimits.m_nuVertexSepMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexSepMax", m_normLimits.m_nuVertexSepMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNHitsMin", m_normLimits.m_parentEndRegionNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNHitsMax", m_normLimits.m_parentEndRegionNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNParticlesMin", m_normLimits.m_parentEndRegionNParticlesMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNParticlesMax", m_normLimits.m_parentEndRegionNParticlesMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionRToWallMin", m_normLimits.m_parentEndRegionRToWallMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionRToWallMax", m_normLimits.m_parentEndRegionRToWallMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexSepMin", m_normLimits.m_vertexSepMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexSepMax", m_normLimits.m_vertexSepMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "DoesChildConnectMin", m_normLimits.m_doesChildConnectMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "DoesChildConnectMax", m_normLimits.m_doesChildConnectMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OvershootDCAMin", m_normLimits.m_overshootDCAMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OvershootDCAMax", m_normLimits.m_overshootDCAMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OvershootLMin", m_normLimits.m_overshootLMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OvershootLMax", m_normLimits.m_overshootLMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChildCPDCAMin", m_normLimits.m_childCPDCAMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChildCPDCAMax", m_normLimits.m_childCPDCAMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ChildCPExtrapDistanceMin", m_normLimits.m_childCPExtrapDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ChildCPExtrapDistanceMax", m_normLimits.m_childCPExtrapDistanceMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ChildCPLRatioMin", m_normLimits.m_childCPLRatioMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ChildCPLRatioMax", m_normLimits.m_childCPLRatioMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitsMin", m_normLimits.m_parentCPNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitsMax", m_normLimits.m_parentCPNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitRatioMin", m_normLimits.m_parentCPNHitRatioMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitRatioMax", m_normLimits.m_parentCPNHitRatioMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPEigenvalueRatioMin", m_normLimits.m_parentCPEigenvalueRatioMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPEigenvalueRatioMax", m_normLimits.m_parentCPEigenvalueRatioMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPOpeningAngleMin", m_normLimits.m_parentCPOpeningAngleMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentCPOpeningAngleMax", m_normLimits.m_parentCPOpeningAngleMax));

    if (!m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackTrackBranchModelName", m_trackTrackBranchModelName));
        m_trackTrackBranchModelName = LArFileHelper::FindFileInPath(m_trackTrackBranchModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackTrackBranchModelName, m_trackTrackBranchModel));

        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackTrackClassifierModelName", m_trackTrackClassifierModelName));
        m_trackTrackClassifierModelName = LArFileHelper::FindFileInPath(m_trackTrackClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackTrackClassifierModelName, m_trackTrackClassifierModel));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackShowerBranchModelName", m_trackShowerBranchModelName));
        m_trackShowerBranchModelName = LArFileHelper::FindFileInPath(m_trackShowerBranchModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackShowerBranchModelName, m_trackShowerBranchModel));

        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackShowerClassifierModelName", m_trackShowerClassifierModelName));
        m_trackShowerClassifierModelName = LArFileHelper::FindFileInPath(m_trackShowerClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackShowerClassifierModelName, m_trackShowerClassifierModel));
    }

    return DLBaseHierarchyTool::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
