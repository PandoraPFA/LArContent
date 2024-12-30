/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.cc
 *
 *  @brief  Implementation of the MLP later tier hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPLaterTierHierarchyTool::MLPLaterTierHierarchyTool() :
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
    m_vertexSepMax(700.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPLaterTierHierarchyTool::Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pNeutrinoPfo, HierarchyPfo &parentHierarchyPfo, 
    HierarchyPfo &childHierarchyPfo)
{
    this->SetDetectorBoundaries();

    if (childHierarchyPfo.GetIsTrack())
    {
        return STATUS_CODE_SUCCESS;
    }
    else
    {
        //////////////////////////////////////////////
        if (this->GetNSpacepoints(parentHierarchyPfo) != 731)
            return STATUS_CODE_SUCCESS;

        if (this->GetNSpacepoints(childHierarchyPfo) != 143)
            return STATUS_CODE_SUCCESS;
        //////////////////////////////////////////////

        // Get shower vertex
        const bool isShowerVertexUpstream(this->IsShowerVertexUpstream(parentHierarchyPfo, childHierarchyPfo));

        // Set network params
        MLPLaterTierNetworkParams upLaterTierNetworkParams, downLaterTierNetworkParams;

        // Get the common params i.e. independent of postulated orientation
        const std::pair<float, float> trackScoreParams(this->GetTrackScoreParams(parentHierarchyPfo, childHierarchyPfo));
        const std::pair<float, float> nSpacepointsParams(this->GetNSpacepointsParams(parentHierarchyPfo, childHierarchyPfo));
        const float separation3D(this->GetSeparation3D(parentHierarchyPfo, childHierarchyPfo));

        // Set common params
        this->SetCommonParams(trackScoreParams, nSpacepointsParams, separation3D, upLaterTierNetworkParams);
        this->SetCommonParams(trackScoreParams, nSpacepointsParams, separation3D, downLaterTierNetworkParams);

        // Get params at parent upstream vertex 
        const StatusCode upStatusCode(this->CalculateNetworkVariables(pAlgorithm, parentHierarchyPfo, childHierarchyPfo, 
            pNeutrinoPfo, true, isShowerVertexUpstream, upLaterTierNetworkParams));

        if (upStatusCode != STATUS_CODE_SUCCESS)
            return upStatusCode;

        // Get params at parent downstream vertex 
        const StatusCode downStatusCode(this->CalculateNetworkVariables(pAlgorithm, parentHierarchyPfo, childHierarchyPfo, 
            pNeutrinoPfo, false, isShowerVertexUpstream, downLaterTierNetworkParams));

        if (downStatusCode != STATUS_CODE_SUCCESS)
            return downStatusCode;

        // // Combine to get primary score
        // const float primaryScore(this->ClassifyTrack(primaryNetworkParamsUp, primaryNetworkParamsDown));

        // hierarchyPfo.SetPrimaryScore(primaryScore);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPLaterTierHierarchyTool::IsShowerVertexUpstream(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo)
{
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
    const HierarchyPfo &childHierarchyPfo)
{
    float parentTrackScore(-999.f), childTrackScore(-999.f);

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
    const HierarchyPfo &childHierarchyPfo)
{
    const float parentNSpacepoints(this->GetNSpacepoints(parentHierarchyPfo));
    const float childNSpacepoints(this->GetNSpacepoints(childHierarchyPfo));

    return std::pair<float, float>({parentNSpacepoints, childNSpacepoints});
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPLaterTierHierarchyTool::GetSeparation3D(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo)
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
    const float separation3D, MLPLaterTierNetworkParams &laterTierNetworkParams)
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
    MLPLaterTierNetworkParams &laterTierNetworkParams)
{
    // Pick out neutrino vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;
 
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    // Pick out the correct pfo vertex/direction
    // Parent POI is a postulate endpoint
    // Child POI is a postulate startpoint
    const CartesianVector parentStart(useUpstreamForParent ? parentHierarchyPfo.GetDownstreamVertex() : parentHierarchyPfo.GetUpstreamVertex());
    const CartesianVector parentEnd(useUpstreamForParent ? parentHierarchyPfo.GetUpstreamVertex() : parentHierarchyPfo.GetDownstreamVertex());
    const CartesianVector parentEndDirection(useUpstreamForParent ? parentHierarchyPfo.GetUpstreamDirection() : parentHierarchyPfo.GetDownstreamDirection());
    const CartesianVector childStart(useUpstreamForChild ? childHierarchyPfo.GetUpstreamVertex() : childHierarchyPfo.GetDownstreamVertex());
    const CartesianVector childStartDirection(useUpstreamForChild ? childHierarchyPfo.GetUpstreamDirection() : childHierarchyPfo.GetDownstreamDirection());

    // Check that we could actually calculate pfo vertex/direction, return if not
    if (!this->IsVectorSet(parentStart)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(parentEnd)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(parentEndDirection)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(childStart)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(childStartDirection)) { return STATUS_CODE_NOT_INITIALIZED; }

    // Set laterTierNetworkParams
    laterTierNetworkParams.m_parentIsPOIClosestToNu = useUpstreamForParent ? 1.f : 0.f;
    laterTierNetworkParams.m_childIsPOIClosestToNu = useUpstreamForChild ? 1.f : 0.f;
    this->SetVertexParams(nuVertex, parentStart, parentEnd, childStart, laterTierNetworkParams);
    this->SetEndRegionParams(pAlgorithm, parentHierarchyPfo.GetPfo(), parentEnd, laterTierNetworkParams);
    this->GetConnectionPoint(parentHierarchyPfo, childHierarchyPfo, parentStart, childStart, childStartDirection);

    /////////////////////////////////////////
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "BEFORE NORMALISATION" << std::endl;
    laterTierNetworkParams.Print();
    std::cout << "-----------------------------------------" << std::endl;
    /////////////////////////////////////////

    // Normalise
    this->NormaliseNetworkParams(laterTierNetworkParams);

    /////////////////////////////////////////
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "AFTER NORMALISATION" << std::endl;
    laterTierNetworkParams.Print();
    std::cout << "-----------------------------------------" << std::endl;
    /////////////////////////////////////////

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::SetVertexParams(const CartesianVector &nuVertex, const CartesianVector &parentStart, 
    const CartesianVector &parentEnd, const CartesianVector &childStart, MLPLaterTierNetworkParams &laterTierNetworkParams)
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
        const float boundaryLow(iAxis == 0 ? m_detectorMinX : iAxis == 1 ? m_detectorMinY : m_detectorMinX);
        const float boundaryHigh(iAxis == 0 ? m_detectorMaxX : iAxis == 1 ? m_detectorMaxY : m_detectorMaxX);

        for (float boundary : {boundaryLow, boundaryHigh})
        {
            const float sep = std::fabs(parentCoord - boundary);
            closestDistance = std::min(closestDistance, sep);
        }
    }

    laterTierNetworkParams.m_parentEndRegionRToWall = closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::GetConnectionPoint(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo, 
    const CartesianVector &parentStart, const CartesianVector &childStart, const CartesianVector &childStartDirection)
{
    CartesianPointVector parentPositions3D;
    LArPfoHelper::GetCoordinateVector(parentHierarchyPfo.GetPfo(), TPC_3D, parentPositions3D);

    const int HALF_WINDOW_LAYERS(20);
    const float pitchU{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)};
    const float slidingFitPitch(std::max({pitchU, pitchV, pitchW}));

    const ThreeDSlidingFitResult slidingFitResult(&parentPositions3D, HALF_WINDOW_LAYERS, slidingFitPitch);

    // Now walk along track, starting from parentStart
    // I want this to be interpreted as l = 0 so have to alter pandora def
    const bool startFromMin((parentStart - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude() >
                            (parentStart - slidingFitResult.GetGlobalMaxLayerPosition()).GetMagnitude());

    const float minL(slidingFitResult.GetLongitudinalDisplacement(slidingFitResult.GetGlobalMinLayerPosition()));
    const float maxL(slidingFitResult.GetLongitudinalDisplacement(slidingFitResult.GetGlobalMaxLayerPosition()));
    const int nSteps(std::floor(std::fabs(maxL - minL) / 1.f));

    float slidingL(startFromMin ? minL : 0.f);

    for (int i = 0; i < nSteps; i++)
    {
        float firstL(slidingL), secondL(slidingL + 1.f);
        float firstEvalL(startFromMin ? firstL : (maxL - firstL));
        float secondEvalL(startFromMin ? secondL : (maxL - secondL));
        CartesianVector firstPosition(0.f, 0.f, 0.f), secondPosition(0.f, 0.f, 0.f);

        if ((slidingFitResult.GetGlobalFitPosition(firstEvalL, firstPosition) != STATUS_CODE_SUCCESS) ||
            (slidingFitResult.GetGlobalFitPosition(secondEvalL, secondPosition) != STATUS_CODE_SUCCESS))
        {
            slidingL += 1.f;
            continue;
        }

        // Extrap child vertex to parent position
        const CartesianVector midpoint((secondPosition + firstPosition) * 0.5f);
        const std::pair<CartesianVector, bool> extrap(this->ExtrapolateChildToParent(midpoint, childStart, childStartDirection));

        std::cout << "######################" << std::endl;
        std::cout << "midpoint: " << midpoint << std::endl;

        if (this->DoesConnect(firstPosition, secondPosition, extrap.first, 50.f))
        {
            std::cout << "HELLO!!" << std::endl;
        }

        slidingL += 1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<CartesianVector, bool> MLPLaterTierHierarchyTool::ExtrapolateChildToParent(const CartesianVector &parentPosition, const CartesianVector &childStart, 
    const CartesianVector &childStartDirection)
{
    const double extrapDistance((parentPosition - childStart).GetDotProduct(childStartDirection));
    const bool backwardsExtrap(extrapDistance < 0.f);
    const CartesianVector extrapPoint(childStart + (childStartDirection * extrapDistance));

    return std::pair<CartesianVector, bool>({extrapPoint, backwardsExtrap});
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPLaterTierHierarchyTool::DoesConnect(const CartesianVector &boundary1, const CartesianVector &boundary2, 
    const CartesianVector &testPoint, const float buffer)
{
    const CartesianVector displacement(boundary2 - boundary1);
    const double segmentLength = displacement.GetMagnitude();
    const CartesianVector unitVector = displacement * (1.f / segmentLength);
    const double l = unitVector.GetDotProduct(testPoint - boundary1);

    std::cout << "l: " << l << std::endl;

    if ((l < 0.f) || (l > segmentLength))
        return false;

    const double t = unitVector.GetCrossProduct(testPoint - boundary1).GetMagnitudeSquared();

    std::cout << "t: " << std::sqrt(t) << std::endl;

    if (t > (buffer * buffer))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPLaterTierHierarchyTool::NormaliseNetworkParams(MLPLaterTierNetworkParams &laterTierNetworkParams)
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPLaterTierHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return MLPBaseHierarchyTool::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
