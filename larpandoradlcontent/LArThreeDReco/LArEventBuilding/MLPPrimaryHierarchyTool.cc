/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.cc
 *
 *  @brief  Implementation of the MLP primary hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"


#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPPrimaryHierarchyTool::MLPPrimaryHierarchyTool() :
    m_detectorMinX(std::numeric_limits<float>::max()),
    m_detectorMaxX(-std::numeric_limits<float>::max()),
    m_detectorMinY(std::numeric_limits<float>::max()),
    m_detectorMaxY(-std::numeric_limits<float>::max()),
    m_detectorMinZ(std::numeric_limits<float>::max()),
    m_detectorMaxZ(-std::numeric_limits<float>::max()),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"}),
    m_vertexRegionRadius(5.f),
    m_extrapolationStepSize(1.f),
    m_nSpacepointsMin(0.f),
    m_nSpacepointsMax(2000.f),
    m_nuSeparationMin(-50.f),
    m_nuSeparationMax(500.f),
    m_vertexRegionNHitsMin(-10.f),
    m_vertexRegionNHitsMax(100.f),
    m_vertexRegionNParticlesMin(-1.f),
    m_vertexRegionNParticlesMax(8.f),
    m_dcaMin(-60.f),
    m_dcaMax(600.f),
    m_connectionExtrapDistanceMin(-700.f),
    m_connectionExtrapDistanceMax(500.f),
    m_parentConnectionDistanceMin(-150.f),
    m_parentConnectionDistanceMax(150.f),
    m_childConnectionDistanceMin(-30.f),
    m_childConnectionDistanceMax(300.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPPrimaryHierarchyTool::Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pNeutrinoPfo, 
    const HierarchyPfoMap &trackPfos, HierarchyPfo &hierarchyPfo)
{
    this->SetDetectorBoundaries();

    if (hierarchyPfo.GetIsTrack())
    {
        // Get params at upstream vertex
        MLPPrimaryNetworkParams primaryNetworkParamsUp;
        primaryNetworkParamsUp.m_isPOIClosestToNu = 1.f;

        const StatusCode statusCodeUp(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, trackPfos, true, primaryNetworkParamsUp));

        if (statusCodeUp != STATUS_CODE_SUCCESS)
            return statusCodeUp;

        // Get params at downstream vertex
        MLPPrimaryNetworkParams primaryNetworkParamsDown;
        primaryNetworkParamsDown.m_isPOIClosestToNu = 0.f;

        const StatusCode statusCodeDown(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, trackPfos, false, primaryNetworkParamsDown));

        if (statusCodeDown != STATUS_CODE_SUCCESS)
            return statusCodeDown;

        // Combine to get primary score
        const float primaryScore(this->ClassifyTrack(primaryNetworkParamsUp, primaryNetworkParamsDown));

        hierarchyPfo.SetPrimaryScore(primaryScore);
    }
    else
    {
        // Get params at upstream vertex 
        MLPPrimaryNetworkParams primaryNetworkParams;
        primaryNetworkParams.m_isPOIClosestToNu = 1.f;

        const StatusCode statusCode(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, trackPfos, true, primaryNetworkParams));

        if (statusCode != STATUS_CODE_SUCCESS)
            return statusCode;

        // Get primary score
        const float primaryScore(this->ClassifyShower(primaryNetworkParams));

        hierarchyPfo.SetPrimaryScore(primaryScore);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetDetectorBoundaries()
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    for (const auto &entry : larTPCMap)
    {
        m_detectorMinX = std::min(m_detectorMinX, entry.second->GetCenterX() - (entry.second->GetWidthX() * 0.5f));
        m_detectorMaxX = std::max(m_detectorMaxX, entry.second->GetCenterX() + (entry.second->GetWidthX() * 0.5f));
        m_detectorMinY = std::min(m_detectorMinY, entry.second->GetCenterY() - (entry.second->GetWidthY() * 0.5f));
        m_detectorMaxY = std::max(m_detectorMaxY, entry.second->GetCenterY() + (entry.second->GetWidthY() * 0.5f));
        m_detectorMinZ = std::min(m_detectorMinZ, entry.second->GetCenterZ() - (entry.second->GetWidthZ() * 0.5f));
        m_detectorMaxZ = std::max(m_detectorMaxZ, entry.second->GetCenterZ() + (entry.second->GetWidthZ() * 0.5f));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPPrimaryHierarchyTool::CalculateNetworkVariables(const Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo, 
    const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const bool useUpstream, MLPPrimaryNetworkParams &primaryNetworkParams)
{
    // Pick out neutrino vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;
 
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    // Pick out the correct pfo vertex/direction
    const CartesianVector particleVertex(useUpstream ? hierarchyPfo.GetUpstreamVertex() : hierarchyPfo.GetDownstreamVertex());
    const CartesianVector particleDirection(useUpstream ? hierarchyPfo.GetUpstreamDirection() : hierarchyPfo.GetDownstreamDirection());

    // Check that we could actually calculate pfo vertex/direction, return if not
    if ((particleVertex.GetZ() < -990.f) || (particleVertex.GetZ() < -990.f))
        return STATUS_CODE_NOT_INITIALIZED;

    if ((particleDirection.GetZ() < -990.f) || (particleDirection.GetZ() < -990.f))
        return STATUS_CODE_NOT_INITIALIZED;

   // Set primaryNetworkParams
    this->SetNSpacepoints(hierarchyPfo.GetPfo(), primaryNetworkParams);
    this->SetNuVertexSep(particleVertex, nuVertex, primaryNetworkParams);
    this->SetVertexRegionParams(pAlgorithm, hierarchyPfo.GetPfo(), particleVertex, primaryNetworkParams);
    this->SetConnectionParams(particleVertex, particleDirection, nuVertex, primaryNetworkParams);
    this->SetContextParams(hierarchyPfo.GetPfo(), particleVertex, particleDirection, nuVertex, trackPfos, primaryNetworkParams);

    // /////////////////////////////////////////
    // std::cout << "-----------------------------------------" << std::endl;
    // std::cout << "BEFORE NORMALISATION" << std::endl;
    // std::cout << "NSpacepoints: " << primaryNetworkParams.m_nSpacepoints << std::endl;
    // std::cout << "NuVertexSep: " << primaryNetworkParams.m_nuSeparation << std::endl;
    // std::cout << "StartRegionNHits: " << primaryNetworkParams.m_vertexRegionNHits << std::endl;
    // std::cout << "StartRegionNParticles: " << primaryNetworkParams.m_vertexRegionNParticles << std::endl;
    // std::cout << "DCA: " << primaryNetworkParams.m_dca << std::endl;
    // std::cout << "ExtrapDistance: " << primaryNetworkParams.m_connectionExtrapDistance << std::endl;
    // std::cout << "parentConnectionDistance: " << primaryNetworkParams.m_parentConnectionDistance << std::endl;
    // std::cout << "childConnectionDistance: " << primaryNetworkParams.m_childConnectionDistance << std::endl;
    // std::cout << "-----------------------------------------" << std::endl;
    // /////////////////////////////////////////

    // Normalise
    this->NormaliseNetworkParams(primaryNetworkParams);

    /////////////////////////////////////////
    // std::cout << "-----------------------------------------" << std::endl;
    // std::cout << "AFTER NORMALISATION" << std::endl;
    // std::cout << "NSpacepoints: " << primaryNetworkParams.m_nSpacepoints << std::endl;
    // std::cout << "NuVertexSep: " << primaryNetworkParams.m_nuSeparation << std::endl;
    // std::cout << "StartRegionNHits: " << primaryNetworkParams.m_vertexRegionNHits << std::endl;
    // std::cout << "StartRegionNParticles: " << primaryNetworkParams.m_vertexRegionNParticles << std::endl;
    // std::cout << "DCA: " << primaryNetworkParams.m_dca << std::endl;
    // std::cout << "ExtrapDistance: " << primaryNetworkParams.m_connectionExtrapDistance << std::endl;
    // std::cout << "parentConnectionDistance: " << primaryNetworkParams.m_parentConnectionDistance << std::endl;
    // std::cout << "childConnectionDistance: " << primaryNetworkParams.m_childConnectionDistance << std::endl;
    // std::cout << "-----------------------------------------" << std::endl;
    /////////////////////////////////////////

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetNSpacepoints(const ParticleFlowObject *const pPfo, MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    primaryNetworkParams.m_nSpacepoints = total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetNuVertexSep(const CartesianVector &particleVertex, const CartesianVector &nuVertex,
    MLPPrimaryNetworkParams &primaryNetworkParams) const
{ 
    const float separation((particleVertex - nuVertex).GetMagnitude());

    primaryNetworkParams.m_nuSeparation = separation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetVertexRegionParams(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, 
    const CartesianVector &particleVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const
{ 
    int hitCount = 0;
    int particleCount = 0;    

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*pAlgorithm, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pOtherPfo : *pPfoList)
        {
            if (pPfo == pOtherPfo)
                continue;

            bool isClose(false);

            CartesianPointVector otherPfoPositions3D;
            LArPfoHelper::GetCoordinateVector(pOtherPfo, TPC_3D, otherPfoPositions3D);

            for (const CartesianVector &otherPfoPosition : otherPfoPositions3D)
            {
                const double sepSq = (otherPfoPosition - particleVertex).GetMagnitudeSquared();

                if (sepSq < (m_vertexRegionRadius * m_vertexRegionRadius))
                {
                    isClose = true;
                    ++hitCount;
                }
            }

            if (isClose)
                ++particleCount;
        }
    }

    primaryNetworkParams.m_vertexRegionNHits = hitCount;
    primaryNetworkParams.m_vertexRegionNParticles = particleCount;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetConnectionParams(const CartesianVector &particleVertex, const CartesianVector &particleDirection, 
    const CartesianVector &nuVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    // Extrapolate particle to nuVertex
    const float extrapDistance((nuVertex - particleVertex).GetDotProduct(particleDirection));
    const CartesianVector extrapolationPoint(particleVertex + (particleDirection * extrapDistance));
    
    primaryNetworkParams.m_dca = (nuVertex - extrapolationPoint).GetMagnitude();
    primaryNetworkParams.m_connectionExtrapDistance = (extrapDistance * (-1.f)); // backwards extrap. should be +ve
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetContextParams(const ParticleFlowObject *const pPfo, const CartesianVector &particleVertex, const CartesianVector &particleDirection, 
    const CartesianVector &nuVertex, const HierarchyPfoMap &trackPfos, MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    // What is our nuVertex separation
    const float nuVertexSepSq((particleVertex - nuVertex).GetMagnitudeSquared());

    bool found(false);
    float parentConnectionDistance(std::numeric_limits<float>::max());
    float childConnectionDistance(std::numeric_limits<float>::max());

    // Look at potential parents
    for (auto& [pTrackPfo, hierarchyTrackPfo] : trackPfos)
    {
        if (pTrackPfo == pPfo)
            continue;

        // Upstream position should be closer to nu vertex
        //const float thisNuVertexSepSq((hierarchyTrackPfo.GetUpstreamVertex() - nuVertex).GetMagnitudeSquared());
        const float thisNuVertexSepSq((hierarchyTrackPfo.GetDownstreamVertex() - nuVertex).GetMagnitudeSquared()); // recreating bug in jupyter code -.-

        if (thisNuVertexSepSq > nuVertexSepSq)
            continue;

        float thisParentConnectionDistance(-999.f), thisChildConnectionDistance(-999.f);

        this->CalculateConnectionDistances(hierarchyTrackPfo.GetDownstreamVertex(), hierarchyTrackPfo.GetDownstreamDirection(), 
            particleVertex, particleDirection, thisParentConnectionDistance, thisChildConnectionDistance);

        if ((thisChildConnectionDistance > 0.f) && (std::fabs(thisParentConnectionDistance) < std::fabs(parentConnectionDistance)))
        {
            found = true;
            parentConnectionDistance = thisParentConnectionDistance;
            childConnectionDistance = thisChildConnectionDistance;
        }
    }

    if (found)
    {
        primaryNetworkParams.m_parentConnectionDistance = parentConnectionDistance;
        primaryNetworkParams.m_childConnectionDistance = childConnectionDistance;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::CalculateConnectionDistances(const CartesianVector &parentEndpoint, const CartesianVector &parentDirection, 
    const CartesianVector &childVertex, const CartesianVector &childDirection, float &parentConnectionDistance, float &childConnectionDistance) const
{
    // Loop things
    float smallestT(std::numeric_limits<float>::max());
    bool isGettingCloser(true);
    bool found(false);
    CartesianVector connectionPoint(-999.f, -999.f, -999.f);

    // start the seed
    CartesianVector extrapolatedPoint(childVertex);

    while (isGettingCloser)
    {
        isGettingCloser = false;
        extrapolatedPoint = extrapolatedPoint + (childDirection * ((-1.f) * m_extrapolationStepSize));

        if (!this->IsInFV(extrapolatedPoint))
            break;

        const float parentT = (parentDirection * (-1.f)).GetCrossProduct((extrapolatedPoint - parentEndpoint)).GetMagnitudeSquared();

        if (parentT < smallestT)
        {
            smallestT = parentT;
            connectionPoint = extrapolatedPoint;
            isGettingCloser = true;
            found = true;
        }
    }        

    // Distance from parent end to connection point
    parentConnectionDistance = found ? (parentDirection * (-1.f)).GetDotProduct(connectionPoint - parentEndpoint) : -999.f; // need to turn parent direction around
    // Distance from child start to connection point
    childConnectionDistance = found ? (childVertex - connectionPoint).GetMagnitude() : -999.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPPrimaryHierarchyTool::IsInFV(const CartesianVector &position) const
{
    if ((position.GetX() < m_detectorMinX) or (position.GetX() > m_detectorMaxX))
        return false;

    if ((position.GetY() < m_detectorMinY) or (position.GetY() > m_detectorMaxY))
        return false;

    if ((position.GetZ() < m_detectorMinZ) or (position.GetZ() > m_detectorMaxZ))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::NormaliseNetworkParams(MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    this->NormaliseNetworkParam(m_nSpacepointsMin, m_nSpacepointsMax, primaryNetworkParams.m_nSpacepoints);
    this->NormaliseNetworkParam(m_nuSeparationMin, m_nuSeparationMax, primaryNetworkParams.m_nuSeparation);
    this->NormaliseNetworkParam(m_vertexRegionNHitsMin, m_vertexRegionNHitsMax, primaryNetworkParams.m_vertexRegionNHits);
    this->NormaliseNetworkParam(m_vertexRegionNParticlesMin, m_vertexRegionNParticlesMax, primaryNetworkParams.m_vertexRegionNParticles);
    this->NormaliseNetworkParam(m_dcaMin, m_dcaMax, primaryNetworkParams.m_dca);
    this->NormaliseNetworkParam(m_connectionExtrapDistanceMin, m_connectionExtrapDistanceMax, primaryNetworkParams.m_connectionExtrapDistance);
    this->NormaliseNetworkParam(m_parentConnectionDistanceMin, m_parentConnectionDistanceMax, primaryNetworkParams.m_parentConnectionDistance);
    this->NormaliseNetworkParam(m_childConnectionDistanceMin, m_childConnectionDistanceMax, primaryNetworkParams.m_childConnectionDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::NormaliseNetworkParam(const float minLimit, const float maxLimit, float &primaryNetworkParam) const
{
    const float interval(std::fabs(minLimit) + std::fabs(maxLimit));

    if (primaryNetworkParam < minLimit)
        primaryNetworkParam = minLimit;

    if (primaryNetworkParam > maxLimit)
        primaryNetworkParam = maxLimit;

    primaryNetworkParam /= interval;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPPrimaryHierarchyTool::ClassifyTrack(const MLPPrimaryNetworkParams &primaryNetworkParamsUp, const MLPPrimaryNetworkParams &primaryNetworkParamsDown)
{
    std::cout << "------------------------" << std::endl;
    std::cout << "Classifying track with " << primaryNetworkParamsUp.m_nSpacepoints << " hits" << std::endl;

    // Invoke branch model for upstream edge
    LArDLHelper::TorchInput inputUp;
    LArDLHelper::InitialiseInput({1, 17}, inputUp);
    inputUp[0][0] = primaryNetworkParamsUp.m_nSpacepoints;
    inputUp[0][1] = primaryNetworkParamsUp.m_nuSeparation;
    inputUp[0][2] = primaryNetworkParamsUp.m_vertexRegionNHits;
    inputUp[0][3] = primaryNetworkParamsUp.m_vertexRegionNParticles;
    inputUp[0][4] = primaryNetworkParamsUp.m_dca;
    inputUp[0][5] = primaryNetworkParamsUp.m_connectionExtrapDistance;
    inputUp[0][6] = primaryNetworkParamsUp.m_isPOIClosestToNu;
    inputUp[0][7] = primaryNetworkParamsUp.m_parentConnectionDistance;
    inputUp[0][8] = primaryNetworkParamsUp.m_childConnectionDistance;
    inputUp[0][9] = primaryNetworkParamsDown.m_nuSeparation;
    inputUp[0][10] = primaryNetworkParamsDown.m_vertexRegionNHits;
    inputUp[0][11] = primaryNetworkParamsDown.m_vertexRegionNParticles;
    inputUp[0][12] = primaryNetworkParamsDown.m_dca;
    inputUp[0][13] = primaryNetworkParamsDown.m_connectionExtrapDistance;
    inputUp[0][14] = primaryNetworkParamsDown.m_isPOIClosestToNu;
    inputUp[0][15] = primaryNetworkParamsDown.m_parentConnectionDistance;
    inputUp[0][16] = primaryNetworkParamsDown.m_childConnectionDistance;
    LArDLHelper::TorchOutput outputUp;
    LArDLHelper::Forward(m_primaryTrackBranchModel, {inputUp}, outputUp);
    auto outputAccessorUp = outputUp.accessor<float, 2>();

    std::cout << "outputUp: " << outputAccessorUp[0][0] << ", " << outputAccessorUp[0][1] << ", " << outputAccessorUp[0][2] << std::endl;

    // Invoke branch model for downstream edge
    LArDLHelper::TorchInput inputDown;
    LArDLHelper::InitialiseInput({1, 17}, inputDown);
    inputDown[0][0] = primaryNetworkParamsDown.m_nSpacepoints;
    inputDown[0][1] = primaryNetworkParamsDown.m_nuSeparation;
    inputDown[0][2] = primaryNetworkParamsDown.m_vertexRegionNHits;
    inputDown[0][3] = primaryNetworkParamsDown.m_vertexRegionNParticles;
    inputDown[0][4] = primaryNetworkParamsDown.m_dca;
    inputDown[0][5] = primaryNetworkParamsDown.m_connectionExtrapDistance;
    inputDown[0][6] = primaryNetworkParamsDown.m_isPOIClosestToNu;
    inputDown[0][7] = primaryNetworkParamsDown.m_parentConnectionDistance;
    inputDown[0][8] = primaryNetworkParamsDown.m_childConnectionDistance;
    inputDown[0][9] = primaryNetworkParamsUp.m_nuSeparation;
    inputDown[0][10] = primaryNetworkParamsUp.m_vertexRegionNHits;
    inputDown[0][11] = primaryNetworkParamsUp.m_vertexRegionNParticles;
    inputDown[0][12] = primaryNetworkParamsUp.m_dca;
    inputDown[0][13] = primaryNetworkParamsUp.m_connectionExtrapDistance;
    inputDown[0][14] = primaryNetworkParamsUp.m_isPOIClosestToNu;
    inputDown[0][15] = primaryNetworkParamsUp.m_parentConnectionDistance;
    inputDown[0][16] = primaryNetworkParamsUp.m_childConnectionDistance;
    LArDLHelper::TorchOutput outputDown;
    LArDLHelper::Forward(m_primaryTrackBranchModel, {inputDown}, outputDown);
    auto outputAccessorDown = outputDown.accessor<float, 2>();

    std::cout << "outputDown: " << outputAccessorDown[0][0] << ", " << outputAccessorDown[0][1] << ", " << outputAccessorDown[0][2] << std::endl;

    // Invoke classifier model
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 6}, input);
    input[0][0] = outputAccessorUp[0][0];
    input[0][1] = outputAccessorUp[0][1];
    input[0][2] = outputAccessorUp[0][2];
    input[0][3] = outputAccessorDown[0][0];
    input[0][4] = outputAccessorDown[0][1];
    input[0][5] = outputAccessorDown[0][2];
    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_primaryTrackClassifierModel, {input}, output);
    auto outputAccessor = output.accessor<float, 2>();

    std::cout << "Track classification score: " << outputAccessor[0][0] << std::endl;
    std::cout << "------------------------" << std::endl;

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPPrimaryHierarchyTool::ClassifyShower(const MLPPrimaryNetworkParams &primaryNetworkParams)
{
    std::cout << "------------------------" << std::endl;
    std::cout << "Classifying shower with " << primaryNetworkParams.m_nSpacepoints << "hits" << std::endl;

    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 9}, input);
    input[0][0] = primaryNetworkParams.m_nSpacepoints;
    input[0][1] = primaryNetworkParams.m_nuSeparation;
    input[0][2] = primaryNetworkParams.m_vertexRegionNHits;
    input[0][3] = primaryNetworkParams.m_vertexRegionNParticles;
    input[0][4] = primaryNetworkParams.m_dca;
    input[0][5] = primaryNetworkParams.m_connectionExtrapDistance;
    input[0][6] = 1.0;
    input[0][7] = primaryNetworkParams.m_parentConnectionDistance;
    input[0][8] = primaryNetworkParams.m_childConnectionDistance;

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_primaryShowerClassifierModel, {input}, output);
    auto outputAccessor = output.accessor<float, 2>();

    std::cout << "Shower classification score: " << outputAccessor[0][0] << std::endl;
    std::cout << "------------------------" << std::endl;
    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPPrimaryHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionRadius", m_vertexRegionRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ExtrapolationStepSize", m_extrapolationStepSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMin", m_nSpacepointsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMax", m_nSpacepointsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuSeparationMin", m_nuSeparationMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuSeparationMax", m_nuSeparationMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNHitsMin", m_vertexRegionNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNHitsMax", m_vertexRegionNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNParticlesMin", m_vertexRegionNParticlesMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNParticlesMax", m_vertexRegionNParticlesMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DCAMin", m_dcaMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DCAMax", m_dcaMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConnectionExtrapDistanceMin", m_connectionExtrapDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConnectionExtrapDistanceMax", m_connectionExtrapDistanceMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ParentConnectionDistanceMin", m_parentConnectionDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ParentConnectionDistanceMax", m_parentConnectionDistanceMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChildConnectionDistanceMin", m_childConnectionDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChildConnectionDistanceMax", m_childConnectionDistanceMax));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryTrackBranchModelName", m_primaryTrackBranchModelName));

    if (!m_primaryTrackBranchModelName.empty())
    {
        m_primaryTrackBranchModelName = LArFileHelper::FindFileInPath(m_primaryTrackBranchModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_primaryTrackBranchModelName, m_primaryTrackBranchModel));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryTrackClassifierModelName", m_primaryTrackClassifierModelName));

    if (!m_primaryTrackClassifierModelName.empty())
    {
        m_primaryTrackClassifierModelName = LArFileHelper::FindFileInPath(m_primaryTrackClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_primaryTrackClassifierModelName, m_primaryTrackClassifierModel));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryShowerClassifierModelName", m_primaryShowerClassifierModelName));

    if (!m_primaryShowerClassifierModelName.empty())
    {
        m_primaryShowerClassifierModelName = LArFileHelper::FindFileInPath(m_primaryShowerClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_primaryShowerClassifierModelName, m_primaryShowerClassifierModel));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
