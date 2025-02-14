/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.cc
 *
 *  @brief  Implementation of the MLP later tier hierarchy tool
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
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPLaterTierHierarchyTool::MLPLaterTierHierarchyTool() :
    MLPBaseHierarchyTool(),
    m_trainingMode(false),
    m_trajectoryStepSize(1.f),
    m_connectionBuffer(50.f),
    m_searchRegion(10.f),
    m_normalise(true),
    m_trackScoreMin(-1.f),
    m_trackScoreMax(1.f),
    m_nSpacepointsMin(0.f),
    m_nSpacepointsMax(2000.f),
    m_separation3DMin(-50.f),
    m_separation3DMax(700.f),
    m_nuVertexSepMin(-100.f),
    m_nuVertexSepMax(750.f),
    m_parentEndRegionNHitsMin(-10.f),
    m_parentEndRegionNHitsMax(80.f),
    m_parentEndRegionNParticlesMin(-1.f),
    m_parentEndRegionNParticlesMax(5.f),
    m_parentEndRegionRToWallMin(-10.f),
    m_parentEndRegionRToWallMax(400.f),
    m_vertexSepMin(-50.f),
    m_vertexSepMax(700.f),
    m_doesChildConnectMin(-1.f),
    m_doesChildConnectMax(1.f),
    m_overshootDCAMin(-700.f),
    m_overshootDCAMax(700.f),
    m_overshootLMin(-100.f),
    m_overshootLMax(700.f),
    m_childCPDCAMin(-5.f),
    m_childCPDCAMax(50.f),
    m_childCPExtrapDistanceMin(-500.f),
    m_childCPExtrapDistanceMax(500.f),
    m_childCPLRatioMin(-5.f),
    m_childCPLRatioMax(30.f),
    m_parentCPNHitsMin(-10.f),
    m_parentCPNHitsMax(100.f),
    m_parentCPNHitRatioMin(-5.f),
    m_parentCPNHitRatioMax(30.f),
    m_parentCPEigenvalueRatioMin(-5.f),
    m_parentCPEigenvalueRatioMax(50.f),
    m_parentCPOpeningAngleMin(-10.f),
    m_parentCPOpeningAngleMax(180.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPLaterTierHierarchyTool::Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &parentHierarchyPfo, 
    const HierarchyPfo &childHierarchyPfo, std::vector<MLPLaterTierNetworkParams> &networkParamVector, float &laterTierScore)
{
    networkParamVector.clear();
    laterTierScore = m_bogusFloat;

    this->SetDetectorBoundaries();

    // Get the common params i.e. independent of postulated orientation
    const std::pair<float, float> trackScoreParams(this->GetTrackScoreParams(parentHierarchyPfo, childHierarchyPfo));
    const std::pair<float, float> nSpacepointsParams(this->GetNSpacepointsParams(parentHierarchyPfo, childHierarchyPfo));
    const float separation3D(this->GetSeparation3D(parentHierarchyPfo, childHierarchyPfo));

    std::vector<bool> childOrientationVector(childHierarchyPfo.GetIsTrack() ? std::vector<bool>({true, false}) : 
        std::vector<bool>({this->IsShowerVertexUpstream(parentHierarchyPfo, childHierarchyPfo)}));

    // Set network params 
    for (const bool &useUpstreamForParent : {true, false})
    {
        for (const bool &useUpstreamForChild : childOrientationVector)
        {
            MLPLaterTierNetworkParams edgeParams;

            // Set common params
            this->SetCommonParams(trackScoreParams, nSpacepointsParams, separation3D, edgeParams);

            // Get orientation dependent params
            const StatusCode statusCode(this->CalculateNetworkVariables(pAlgorithm, parentHierarchyPfo, childHierarchyPfo, 
                pNeutrinoPfo, useUpstreamForParent, useUpstreamForChild, edgeParams));

            if (statusCode != STATUS_CODE_SUCCESS)
                return statusCode;

            // Add these to our output vector
            networkParamVector.emplace_back(edgeParams);
        }
    }

    // Now run the model!
    if (!m_trainingMode)
    {
        if (childHierarchyPfo.GetIsTrack())
            laterTierScore = this->ClassifyTrackTrack(networkParamVector.at(0), networkParamVector.at(1), networkParamVector.at(2), networkParamVector.at(3));
        else
            laterTierScore = this->ClassifyTrackShower(networkParamVector.at(0), networkParamVector.at(1));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPLaterTierHierarchyTool::IsShowerVertexUpstream(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const
{
    const bool isUpstreamSet(this->IsVectorSet(childHierarchyPfo.GetUpstreamVertex()));
    const bool isDownstreamSet(this->IsVectorSet(childHierarchyPfo.GetDownstreamVertex()));

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
        const float thisUpstreamSepSq((parentPosition3D - childHierarchyPfo.GetUpstreamVertex()).GetMagnitudeSquared());
        const float thisDownstreamSepSq((parentPosition3D - childHierarchyPfo.GetDownstreamVertex()).GetMagnitudeSquared());

        upstreamSepSq = std::min(upstreamSepSq, thisUpstreamSepSq);
        downstreamSepSq = std::min(downstreamSepSq, thisDownstreamSepSq);
    }

    return (upstreamSepSq < downstreamSepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, float> MLPLaterTierHierarchyTool::GetTrackScoreParams(const HierarchyPfo &parentHierarchyPfo, 
    const HierarchyPfo &childHierarchyPfo) const
{
    float parentTrackScore(m_bogusFloat), childTrackScore(m_bogusFloat);

    const PropertiesMap &parentMeta(parentHierarchyPfo.GetPfo()->GetPropertiesMap());
    const PropertiesMap &childMeta(childHierarchyPfo.GetPfo()->GetPropertiesMap());

    if (parentMeta.find("TrackScore") != parentMeta.end())
        parentTrackScore = parentMeta.at("TrackScore");

    if (childMeta.find("TrackScore") != childMeta.end())
        childTrackScore = childMeta.at("TrackScore");

    return std::pair<float, float>({parentTrackScore, childTrackScore});
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, float> MLPLaterTierHierarchyTool::GetNSpacepointsParams(const HierarchyPfo &parentHierarchyPfo, 
    const HierarchyPfo &childHierarchyPfo) const
{
    const float parentNSpacepoints(this->GetNSpacepoints(parentHierarchyPfo));
    const float childNSpacepoints(this->GetNSpacepoints(childHierarchyPfo));

    return std::pair<float, float>({parentNSpacepoints, childNSpacepoints});
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPLaterTierHierarchyTool::GetSeparation3D(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const
{
    CartesianPointVector parentPositions3D, childPositions3D;
    LArPfoHelper::GetCoordinateVector(parentHierarchyPfo.GetPfo(), TPC_3D, parentPositions3D);
    LArPfoHelper::GetCoordinateVector(childHierarchyPfo.GetPfo(), TPC_3D, childPositions3D);

    float sepSq(std::numeric_limits<float>::max());

    for (const CartesianVector &parentPosition3D : parentPositions3D)
    {
        for (const CartesianVector &childPosition3D : childPositions3D)
        {
            const float thisSepSq((parentPosition3D - childPosition3D).GetMagnitudeSquared());

            sepSq = std::min(thisSepSq, sepSq);
        }
    }

    return std::sqrt(sepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetCommonParams(const std::pair<float, float> &trackScoreParams, const std::pair<float, float> &nSpacepointsParams, 
    const float separation3D, MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    laterTierNetworkParams.m_parentTrackScore = trackScoreParams.first;
    laterTierNetworkParams.m_childTrackScore = trackScoreParams.second;
    laterTierNetworkParams.m_parentNSpacepoints = nSpacepointsParams.first;
    laterTierNetworkParams.m_childNSpacepoints = nSpacepointsParams.second;
    laterTierNetworkParams.m_separation3D = separation3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPLaterTierHierarchyTool::CalculateNetworkVariables(const Algorithm *const pAlgorithm, const HierarchyPfo &parentHierarchyPfo, 
    const HierarchyPfo &childHierarchyPfo, const ParticleFlowObject *const pNeutrinoPfo, const bool useUpstreamForParent, const bool useUpstreamForChild, 
    MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    // Pick out neutrino vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;
 
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    // Pick out the correct pfo vertex/direction
    // Parent POI is a postulate endpoint
    // Child POI is a postulate startpoint
    const CartesianVector &parentStart(useUpstreamForParent ? parentHierarchyPfo.GetDownstreamVertex() : parentHierarchyPfo.GetUpstreamVertex());
    const CartesianVector &parentStartDirection(useUpstreamForParent ? parentHierarchyPfo.GetDownstreamDirection() : parentHierarchyPfo.GetUpstreamDirection());
    const CartesianVector &parentEnd(useUpstreamForParent ? parentHierarchyPfo.GetUpstreamVertex() : parentHierarchyPfo.GetDownstreamVertex());
    const CartesianVector &parentEndDirection(useUpstreamForParent ? parentHierarchyPfo.GetUpstreamDirection() : parentHierarchyPfo.GetDownstreamDirection());
    const CartesianVector &childStart(useUpstreamForChild ? childHierarchyPfo.GetUpstreamVertex() : childHierarchyPfo.GetDownstreamVertex());
    const CartesianVector &childStartDirection(useUpstreamForChild ? childHierarchyPfo.GetUpstreamDirection() : childHierarchyPfo.GetDownstreamDirection());

    // Check that we could actually calculate pfo vertex/direction, return if not
    if (!this->IsVectorSet(parentStart)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(parentStartDirection)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(parentEnd)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(parentEndDirection)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(childStart)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(childStartDirection)) { return STATUS_CODE_NOT_INITIALIZED; }

    // Set laterTierNetworkParams
    laterTierNetworkParams.m_parentIsPOIClosestToNu = useUpstreamForParent ? 1.f : 0.f;
    laterTierNetworkParams.m_childIsPOIClosestToNu = useUpstreamForChild ? 1.f : 0.f;
    this->SetVertexParams(nuVertex, parentStart, parentEnd, childStart, laterTierNetworkParams);
    this->SetEndRegionParams(pAlgorithm, parentHierarchyPfo.GetPfo(), parentEnd, laterTierNetworkParams);
    this->SetConnectionParams(parentHierarchyPfo, parentStart, childStart, childStartDirection, laterTierNetworkParams);
    this->SetOvershootParams(parentStart, parentStartDirection, parentEnd, parentEndDirection, childStart, childStartDirection, laterTierNetworkParams);
    this->SetParentConnectionPointVars(parentHierarchyPfo, laterTierNetworkParams);

    /////////////////////////////////////////
    // std::cout << "-----------------------------------------" << std::endl;
    // std::cout << "BEFORE NORMALISATION" << std::endl;
    // laterTierNetworkParams.Print();
    // std::cout << "-----------------------------------------" << std::endl;
    /////////////////////////////////////////

    // Normalise
    if (m_normalise)
        this->NormaliseNetworkParams(laterTierNetworkParams);

    /////////////////////////////////////////
    // std::cout << "-----------------------------------------" << std::endl;
    // std::cout << "AFTER NORMALISATION" << std::endl;
    // laterTierNetworkParams.Print();
    // std::cout << "-----------------------------------------" << std::endl;
    /////////////////////////////////////////

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetVertexParams(const CartesianVector &nuVertex, const CartesianVector &parentStart, 
    const CartesianVector &parentEnd, const CartesianVector &childStart, MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    laterTierNetworkParams.m_parentNuVertexSep = ((parentStart - nuVertex).GetMagnitude());
    laterTierNetworkParams.m_childNuVertexSep = ((childStart - nuVertex).GetMagnitude());
    laterTierNetworkParams.m_vertexSeparation = ((parentEnd - childStart).GetMagnitude());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetEndRegionParams(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pParentPfo, 
    const CartesianVector &parentEnd, MLPLaterTierNetworkParams &laterTierNetworkParams) const
{ 
    std::pair<float, float> endRegionParams(this->GetParticleInfoAboutPfoPosition(pAlgorithm, pParentPfo, parentEnd));

    laterTierNetworkParams.m_parentEndRegionNHits = endRegionParams.first;
    laterTierNetworkParams.m_parentEndRegionNParticles = endRegionParams.second;

    this->SetEndRegionRToWall(parentEnd, laterTierNetworkParams);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetEndRegionRToWall(const CartesianVector &parentEnd, MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    float closestDistance = std::numeric_limits<float>::max();

    for (int iAxis : {0, 1, 2})
    {
        const float parentCoord(iAxis == 0 ? parentEnd.GetX() : iAxis == 1 ? parentEnd.GetY() : parentEnd.GetZ());
        const float boundaryLow(iAxis == 0 ? m_detectorMinX : iAxis == 1 ? m_detectorMinY : m_detectorMinZ);
        const float boundaryHigh(iAxis == 0 ? m_detectorMaxX : iAxis == 1 ? m_detectorMaxY : m_detectorMaxZ);

        for (float boundary : {boundaryLow, boundaryHigh})
        {
            const float sep = std::fabs(parentCoord - boundary);
            closestDistance = std::min(closestDistance, sep);
        }
    }

    laterTierNetworkParams.m_parentEndRegionRToWall = closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetConnectionParams(const HierarchyPfo &parentHierarchyPfo, 
    const CartesianVector &parentStart, const CartesianVector &childStart, const CartesianVector &childStartDirection, 
    MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    // Get fit of parent
    const ThreeDSlidingFitResult &slidingFitResult(parentHierarchyPfo.GetSlidingFitResult());

    // Now walk along track, starting from parentStart
    // I want this to be interpreted as l = 0 so have to alter pandora def
    const bool startFromMin((parentStart - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitudeSquared() <
                            (parentStart - slidingFitResult.GetGlobalMaxLayerPosition()).GetMagnitudeSquared());
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
        CartesianVector firstPosition(m_bogusFloat, m_bogusFloat, m_bogusFloat);
        CartesianVector secondPosition(m_bogusFloat, m_bogusFloat, m_bogusFloat);

        if ((slidingFitResult.GetGlobalFitPosition(firstEvalL, firstPosition) != STATUS_CODE_SUCCESS) ||
            (slidingFitResult.GetGlobalFitPosition(secondEvalL, secondPosition) != STATUS_CODE_SUCCESS))
        {
            slidingL += m_trajectoryStepSize;
            continue;
        }

        // Get trajectory position
        const CartesianVector midPoint((secondPosition + firstPosition) * 0.5f);

        // Extrap child vertex to parent position
        const std::pair<CartesianVector, bool> extrap(this->ExtrapolateChildToParent(midPoint, childStart, childStartDirection));
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
                laterTierNetworkParams.m_childCPExtrapDistance = (extrap.first - childStart).GetMagnitude() * (extrap.second ? 1.f : (-1.f));
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

std::pair<CartesianVector, bool> MLPLaterTierHierarchyTool::ExtrapolateChildToParent(const CartesianVector &parentPosition, const CartesianVector &childStart, 
    const CartesianVector &childStartDirection) const
{
    const double extrapDistance((parentPosition - childStart).GetDotProduct(childStartDirection));
    const bool backwardsExtrap(extrapDistance < 0.f);
    const CartesianVector extrapPoint(childStart + (childStartDirection * extrapDistance));

    return std::pair<CartesianVector, bool>({extrapPoint, backwardsExtrap});
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPLaterTierHierarchyTool::DoesConnect(const CartesianVector &boundary1, const CartesianVector &boundary2, 
    const CartesianVector &testPoint) const
{
    const CartesianVector displacement(boundary2 - boundary1);
    const double segmentLength = displacement.GetMagnitude();
    const CartesianVector unitVector = displacement * (1.f / segmentLength);
    const double l = unitVector.GetDotProduct(testPoint - boundary1);

    if ((l < 0.f) || (l > segmentLength))
        return false;

    const double t = unitVector.GetCrossProduct(testPoint - boundary1).GetMagnitudeSquared();

    if (t > (m_connectionBuffer * m_connectionBuffer))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetOvershootParams(const CartesianVector &parentStart, const CartesianVector &parentStartDirection, 
    const CartesianVector &parentEnd, const CartesianVector &parentEndDirection, const CartesianVector &childStart, 
    const CartesianVector &childStartDirection, MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    if (std::fabs(laterTierNetworkParams.m_doesChildConnect - 1.f) < std::numeric_limits<float>::epsilon())
        return;

    // Let's start at the very beginning, a very good place to start
    const std::pair<CartesianVector, bool> extrapStart(this->ExtrapolateChildToParent(parentStart, childStart, childStartDirection));
    laterTierNetworkParams.m_overshootStartDCA = (extrapStart.first - parentStart).GetMagnitude() * (extrapStart.second ? 1.f : (-1.f));
    laterTierNetworkParams.m_overshootStartL = std::fabs((childStart - parentStart).GetDotProduct(parentStartDirection));

    // Now end
    const std::pair<CartesianVector, bool> extrapEnd(this->ExtrapolateChildToParent(parentEnd, childStart, childStartDirection));
    laterTierNetworkParams.m_overshootEndDCA = (extrapEnd.first - parentEnd).GetMagnitude() * (extrapEnd.second ? 1.f : (-1.f));
    laterTierNetworkParams.m_overshootEndL = std::fabs((childStart - parentEnd).GetDotProduct(parentEndDirection));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetParentConnectionPointVars(const HierarchyPfo &parentHierarchyPfo, 
    MLPLaterTierNetworkParams &laterTierNetworkParams) const
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
        const float thisL = laterTierNetworkParams.m_connectionDirection.GetDotProduct(displacement);
        CartesianPointVector &group(thisL > 0.f ? downstreamGroup : upstreamGroup);
        group.emplace_back(parentPosition);
    }

    laterTierNetworkParams.m_parentCPNUpstreamHits = upstreamGroup.size();
    laterTierNetworkParams.m_parentCPNDownstreamHits = downstreamGroup.size();

    if (upstreamGroup.empty() || downstreamGroup.empty())
        return;

    laterTierNetworkParams.m_parentCPNHitRatio = 
        laterTierNetworkParams.m_parentCPNDownstreamHits / laterTierNetworkParams.m_parentCPNUpstreamHits;

    // Now PCA magic
    try
    {
        CartesianVector centroidUp(m_bogusFloat, m_bogusFloat, m_bogusFloat);
        CartesianVector centroidDown(m_bogusFloat, m_bogusFloat, m_bogusFloat);
        LArPcaHelper::EigenVectors eigenVecsUp, eigenVecsDown;
        LArPcaHelper::EigenValues eigenValuesUp(m_bogusFloat, m_bogusFloat, m_bogusFloat);
        LArPcaHelper::EigenValues  eigenValuesDown(m_bogusFloat, m_bogusFloat, m_bogusFloat);

        LArPcaHelper::RunPca(upstreamGroup, centroidUp, eigenValuesUp, eigenVecsUp);
        LArPcaHelper::RunPca(downstreamGroup, centroidDown, eigenValuesDown, eigenVecsDown);

        // Get opening angle from first eigenvectors (this is the longitudinal one) - straight would be around 180
        laterTierNetworkParams.m_parentCPOpeningAngle = eigenVecsUp.at(0).GetOpeningAngle(eigenVecsDown.at(0)) * (180.f / 3.14);

        // Get average transverse eigenvalues, get ratio
        const float avTransEVaUp((eigenValuesUp.GetY() + eigenValuesUp.GetZ()) * 0.5f);
        const float avTransEVaDown((eigenValuesDown.GetY() + eigenValuesDown.GetZ()) * 0.5f);

        laterTierNetworkParams.m_parentCPEigenvalueRatio = avTransEVaUp < std::numeric_limits<double>::epsilon() ? 0.f : 
            avTransEVaDown / avTransEVaUp;
    }
    catch (...) { return; }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::NormaliseNetworkParams(MLPLaterTierNetworkParams &laterTierNetworkParams) const
{
    this->NormaliseNetworkParam(m_trackScoreMin, m_trackScoreMax, laterTierNetworkParams.m_parentTrackScore);
    this->NormaliseNetworkParam(m_trackScoreMin, m_trackScoreMax, laterTierNetworkParams.m_childTrackScore);
    this->NormaliseNetworkParam(m_nSpacepointsMin, m_nSpacepointsMax, laterTierNetworkParams.m_parentNSpacepoints);
    this->NormaliseNetworkParam(m_nSpacepointsMin, m_nSpacepointsMax, laterTierNetworkParams.m_childNSpacepoints);
    this->NormaliseNetworkParam(m_separation3DMin, m_separation3DMax, laterTierNetworkParams.m_separation3D);
    this->NormaliseNetworkParam(m_nuVertexSepMin, m_nuVertexSepMax, laterTierNetworkParams.m_parentNuVertexSep); 
    this->NormaliseNetworkParam(m_nuVertexSepMin, m_nuVertexSepMax, laterTierNetworkParams.m_childNuVertexSep);
    this->NormaliseNetworkParam(m_parentEndRegionNHitsMin, m_parentEndRegionNHitsMax, laterTierNetworkParams.m_parentEndRegionNHits);
    this->NormaliseNetworkParam(m_parentEndRegionNParticlesMin, m_parentEndRegionNParticlesMax, laterTierNetworkParams.m_parentEndRegionNParticles);
    this->NormaliseNetworkParam(m_parentEndRegionRToWallMin, m_parentEndRegionRToWallMax, laterTierNetworkParams.m_parentEndRegionRToWall);
    this->NormaliseNetworkParam(m_vertexSepMin, m_vertexSepMax, laterTierNetworkParams.m_vertexSeparation);
    this->NormaliseNetworkParam(m_doesChildConnectMin, m_doesChildConnectMax, laterTierNetworkParams.m_doesChildConnect);
    this->NormaliseNetworkParam(m_overshootDCAMin, m_overshootDCAMax, laterTierNetworkParams.m_overshootStartDCA);
    this->NormaliseNetworkParam(m_overshootLMin, m_overshootLMax, laterTierNetworkParams.m_overshootStartL);
    this->NormaliseNetworkParam(m_overshootDCAMin, m_overshootDCAMax, laterTierNetworkParams.m_overshootEndDCA);
    this->NormaliseNetworkParam(m_overshootLMin, m_overshootLMax, laterTierNetworkParams.m_overshootEndL);
    this->NormaliseNetworkParam(m_childCPDCAMin, m_childCPDCAMax, laterTierNetworkParams.m_childCPDCA);
    this->NormaliseNetworkParam(m_childCPExtrapDistanceMin, m_childCPExtrapDistanceMax, laterTierNetworkParams.m_childCPExtrapDistance);
    this->NormaliseNetworkParam(m_childCPLRatioMin, m_childCPLRatioMax, laterTierNetworkParams.m_childCPLRatio);
    this->NormaliseNetworkParam(m_parentCPNHitsMin, m_parentCPNHitsMax, laterTierNetworkParams.m_parentCPNUpstreamHits);
    this->NormaliseNetworkParam(m_parentCPNHitsMin, m_parentCPNHitsMax, laterTierNetworkParams.m_parentCPNDownstreamHits);
    this->NormaliseNetworkParam(m_parentCPNHitRatioMin, m_parentCPNHitRatioMax, laterTierNetworkParams.m_parentCPNHitRatio);
    this->NormaliseNetworkParam(m_parentCPEigenvalueRatioMin, m_parentCPEigenvalueRatioMax, laterTierNetworkParams.m_parentCPEigenvalueRatio);
    this->NormaliseNetworkParam(m_parentCPOpeningAngleMin, m_parentCPOpeningAngleMax, laterTierNetworkParams.m_parentCPOpeningAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPLaterTierHierarchyTool::ClassifyTrackTrack(const MLPLaterTierNetworkParams &edgeParamsUpUp, const MLPLaterTierNetworkParams &edgeParamsUpDown, 
    const MLPLaterTierNetworkParams &edgeParamsDownUp, const MLPLaterTierNetworkParams &edgeParamsDownDown)
{
    ////////////////////////////////////////////////////////////
    // std::cout << "------------------------" << std::endl;
    // std::cout << "Classifying track-track edge with " << edgeParamsUpUp.m_parentNSpacepoints << " parent spacepoints and " << 
    //     edgeParamsUpUp.m_childNSpacepoints << " child spacepoints" << std::endl;
    ////////////////////////////////////////////////////////////

    // Invoke branch model for each edge
    const FloatVector outputUpUp(this->ClassifyTrackTrackEdge(edgeParamsUpUp, edgeParamsUpDown, edgeParamsDownUp, edgeParamsDownDown));
    const FloatVector outputUpDown(this->ClassifyTrackTrackEdge(edgeParamsUpDown, edgeParamsDownUp, edgeParamsDownDown, edgeParamsUpUp));
    const FloatVector outputDownUp(this->ClassifyTrackTrackEdge(edgeParamsDownUp, edgeParamsDownDown, edgeParamsUpUp, edgeParamsUpDown));
    const FloatVector outputDownDown(this->ClassifyTrackTrackEdge(edgeParamsDownDown, edgeParamsUpUp, edgeParamsUpDown, edgeParamsDownUp));

    ////////////////////////////////////////////////////////////
    // std::cout << "outputUpUp: " << outputUpUp.at(0) << ", " << outputUpUp.at(1) << ", " << outputUpUp.at(2) << std::endl;
    // std::cout << "outputUpDown: " << outputUpDown.at(0) << ", " << outputUpDown.at(1) << ", " << outputUpDown.at(2) << std::endl;
    // std::cout << "outputDownUp: " << outputDownUp.at(0) << ", " << outputDownUp.at(1) << ", " << outputDownUp.at(2) << std::endl;
    // std::cout << "outputDownDown: " << outputDownDown.at(0) << ", " << outputDownDown.at(1) << ", " << outputDownDown.at(2) << std::endl;
    ////////////////////////////////////////////////////////////

    // Invoke classifier model for final output
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 12}, input);

    int insertIndex = 0;

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

    ////////////////////////////////////////////////////////////
    // std::cout << "Track-track classification score: " << outputAccessor[0][0] << std::endl;
    // std::cout << "------------------------" << std::endl;
    ////////////////////////////////////////////////////////////

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector MLPLaterTierHierarchyTool::ClassifyTrackTrackEdge(const MLPLaterTierNetworkParams &edgeParams, 
    const MLPLaterTierNetworkParams &otherEdgeParams1, const MLPLaterTierNetworkParams &otherEdgeParams2, const MLPLaterTierNetworkParams &otherEdgeParams3)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 89}, input);

    int insertIndex(this->AddToInput(0, edgeParams.GetCommonParamsForModel(), input));
    insertIndex = this->AddToInput(insertIndex, edgeParams.GetOrientationParamsForModel(), input);
    insertIndex = this->AddToInput(insertIndex, otherEdgeParams1.GetOrientationParamsForModel(), input);
    insertIndex = this->AddToInput(insertIndex, otherEdgeParams2.GetOrientationParamsForModel(), input);
    insertIndex = this->AddToInput(insertIndex, otherEdgeParams3.GetOrientationParamsForModel(), input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_trackTrackBranchModel, {input}, output);

    const torch::TensorAccessor<float, 2> outputAccessor(output.accessor<float, 2>());

    return {outputAccessor[0][0], outputAccessor[0][1], outputAccessor[0][2]};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPLaterTierHierarchyTool::ClassifyTrackShower(const MLPLaterTierNetworkParams &edgeParamsUp, const MLPLaterTierNetworkParams &edgeParamsDown)
{
    ////////////////////////////////////////////////////////////
    // std::cout << "------------------------" << std::endl;
    // std::cout << "Classifying track-shower edge with " << edgeParamsUp.m_parentNSpacepoints << " parent spacepoints and " << 
    //     edgeParamsUp.m_childNSpacepoints << " child spacepoints" << std::endl;
    ////////////////////////////////////////////////////////////

    // Invoke branch model for each edge
    const FloatVector outputUp(this->ClassifyTrackShowerEdge(edgeParamsUp, edgeParamsDown));
    const FloatVector outputDown(this->ClassifyTrackShowerEdge(edgeParamsDown, edgeParamsUp));

    ////////////////////////////////////////////////////////////
    // std::cout << "outputUp: " << outputUp.at(0) << ", " << outputUp.at(1) << ", " << outputUp.at(2) << std::endl;
    // std::cout << "outputDown: " << outputDown.at(0) << ", " << outputDown.at(1) << ", " << outputDown.at(2) << std::endl;
    ////////////////////////////////////////////////////////////

    // Invoke classifier model for final output
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 6}, input);

    int insertIndex = 0;

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

    ////////////////////////////////////////////////////////////
    // std::cout << "Track-shower classification score: " << outputAccessor[0][0] << std::endl;
    // std::cout << "------------------------" << std::endl;
    ////////////////////////////////////////////////////////////

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector MLPLaterTierHierarchyTool::ClassifyTrackShowerEdge(const MLPLaterTierNetworkParams &edgeParams, 
    const MLPLaterTierNetworkParams &otherEdgeParams)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 47}, input);

    int insertIndex(this->AddToInput(0, edgeParams.GetCommonParamsForModel(), input));
    insertIndex = this->AddToInput(insertIndex, edgeParams.GetOrientationParamsForModel(), input);
    insertIndex = this->AddToInput(insertIndex, otherEdgeParams.GetOrientationParamsForModel(), input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_trackShowerBranchModel, {input}, output);

    const torch::TensorAccessor<float, 2> outputAccessor(output.accessor<float, 2>());

    return {outputAccessor[0][0], outputAccessor[0][1], outputAccessor[0][2]};
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPLaterTierHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrajectoryStepSize", m_trajectoryStepSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConnectionBuffer", m_connectionBuffer));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SearchRegion", m_searchRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Normalise", m_normalise));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TrackScoreMin", m_trackScoreMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TrackScoreMax", m_trackScoreMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "NSpacepointsMin", m_nSpacepointsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "NSpacepointsMax", m_nSpacepointsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "Separation3DMin", m_separation3DMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "Separation3DMax", m_separation3DMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "NuVertexSepMin", m_nuVertexSepMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "NuVertexSepMax", m_nuVertexSepMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNHitsMin", m_parentEndRegionNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNHitsMax", m_parentEndRegionNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNParticlesMin", m_parentEndRegionNParticlesMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionNParticlesMax", m_parentEndRegionNParticlesMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionRToWallMin", m_parentEndRegionRToWallMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentEndRegionRToWallMax", m_parentEndRegionRToWallMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "VertexSepMin", m_vertexSepMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "VertexSepMax", m_vertexSepMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "DoesChildConnectMin", m_doesChildConnectMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "DoesChildConnectMax", m_doesChildConnectMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "OvershootDCAMin", m_overshootDCAMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "OvershootDCAMax", m_overshootDCAMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "OvershootLMin", m_overshootLMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "OvershootLMax", m_overshootLMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ChildCPDCAMin", m_childCPDCAMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ChildCPDCAMax", m_childCPDCAMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ChildCPExtrapDistanceMin", m_childCPExtrapDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ChildCPExtrapDistanceMax", m_childCPExtrapDistanceMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ChildCPLRatioMin", m_childCPLRatioMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ChildCPLRatioMax", m_childCPLRatioMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitsMin", m_parentCPNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitsMax", m_parentCPNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitRatioMin", m_parentCPNHitRatioMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPNHitRatioMax", m_parentCPNHitRatioMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPEigenvalueRatioMin", m_parentCPEigenvalueRatioMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPEigenvalueRatioMax", m_parentCPEigenvalueRatioMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPOpeningAngleMin", m_parentCPOpeningAngleMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ParentCPOpeningAngleMax", m_parentCPOpeningAngleMax));

    if (!m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackTrackBranchModelName", m_trackTrackBranchModelName));        
        m_trackTrackBranchModelName = LArFileHelper::FindFileInPath(m_trackTrackBranchModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackTrackBranchModelName, m_trackTrackBranchModel));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackTrackClassifierModelName", m_trackTrackClassifierModelName));        
        m_trackTrackClassifierModelName = LArFileHelper::FindFileInPath(m_trackTrackClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackTrackClassifierModelName, m_trackTrackClassifierModel));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackShowerBranchModelName", m_trackShowerBranchModelName));
        m_trackShowerBranchModelName = LArFileHelper::FindFileInPath(m_trackShowerBranchModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackShowerBranchModelName, m_trackShowerBranchModel));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackShowerClassifierModelName", m_trackShowerClassifierModelName));        
        m_trackShowerClassifierModelName = LArFileHelper::FindFileInPath(m_trackShowerClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_trackShowerClassifierModelName, m_trackShowerClassifierModel));
    }

    return MLPBaseHierarchyTool::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
