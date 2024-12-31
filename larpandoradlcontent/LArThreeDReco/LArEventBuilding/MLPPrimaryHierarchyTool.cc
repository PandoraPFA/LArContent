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
    MLPBaseHierarchyTool(),
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
        // Set network params
        MLPPrimaryNetworkParams primaryNetworkParamsUp, primaryNetworkParamsDown;

        const StatusCode statusCodeUp(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, 
            trackPfos, true, primaryNetworkParamsUp));

        if (statusCodeUp != STATUS_CODE_SUCCESS)
            return statusCodeUp;

        const StatusCode statusCodeDown(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, 
            trackPfos, false, primaryNetworkParamsDown));

        if (statusCodeDown != STATUS_CODE_SUCCESS)
            return statusCodeDown;

        // Now run the model!
        const float primaryScore(this->ClassifyTrack(primaryNetworkParamsUp, primaryNetworkParamsDown));

        hierarchyPfo.SetPrimaryScore(primaryScore);
    }
    else
    {
        // Set network params
        MLPPrimaryNetworkParams primaryNetworkParams;

        const StatusCode statusCode(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, 
            trackPfos, true, primaryNetworkParams));

        if (statusCode != STATUS_CODE_SUCCESS)
            return statusCode;

        // Now run the model!
        const float primaryScore(this->ClassifyShower(primaryNetworkParams));

        hierarchyPfo.SetPrimaryScore(primaryScore);
    }

    return STATUS_CODE_SUCCESS;
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
    if (!this->IsVectorSet(particleVertex)) { return STATUS_CODE_NOT_INITIALIZED; }
    if (!this->IsVectorSet(particleDirection)) { return STATUS_CODE_NOT_INITIALIZED; }

    // Set primaryNetworkParams
    primaryNetworkParams.m_isPOIClosestToNu = useUpstream ? 1.f : 0.f;
    primaryNetworkParams.m_nSpacepoints = this->GetNSpacepoints(hierarchyPfo);
    primaryNetworkParams.m_nuSeparation = (particleVertex - nuVertex).GetMagnitude();
    this->SetVertexRegionParams(pAlgorithm, hierarchyPfo.GetPfo(), particleVertex, primaryNetworkParams);
    this->SetConnectionParams(particleVertex, particleDirection, nuVertex, primaryNetworkParams);
    this->SetContextParams(hierarchyPfo.GetPfo(), particleVertex, particleDirection, nuVertex, trackPfos, primaryNetworkParams);

    // /////////////////////////////////////////
    // std::cout << "-----------------------------------------" << std::endl;
    // std::cout << "BEFORE NORMALISATION" << std::endl;
    // primaryNetworkParams.Print();
    // std::cout << "-----------------------------------------" << std::endl;
    // /////////////////////////////////////////

    // Normalise
    this->NormaliseNetworkParams(primaryNetworkParams);

    /////////////////////////////////////////
    // std::cout << "-----------------------------------------" << std::endl;
    // std::cout << "AFTER NORMALISATION" << std::endl;
    // primaryNetworkParams.Print();
    // std::cout << "-----------------------------------------" << std::endl;
    /////////////////////////////////////////

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetVertexRegionParams(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, 
    const CartesianVector &particleVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const
{ 
    std::pair<float, float> vertexRegionParams(this->GetParticleInfoAboutPfoPosition(pAlgorithm, pPfo, particleVertex));

    primaryNetworkParams.m_vertexRegionNHits = vertexRegionParams.first;
    primaryNetworkParams.m_vertexRegionNParticles = vertexRegionParams.second;
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

float MLPPrimaryHierarchyTool::ClassifyTrack(const MLPPrimaryNetworkParams &edgeParamsUp, const MLPPrimaryNetworkParams &edgeParamsDown)
{
    ////////////////////////////////////////////////////////////
    std::cout << "------------------------" << std::endl;
    std::cout << "Classifying track with " << edgeParamsUp.m_nSpacepoints << " hits" << std::endl;
    ////////////////////////////////////////////////////////////

    // Invoke branch model for each edge
    torch::TensorAccessor<float, 2> outputAccessorUp(this->ClassifyTrackEdge(edgeParamsUp, edgeParamsDown));
    torch::TensorAccessor<float, 2> outputAccessorDown(this->ClassifyTrackEdge(edgeParamsDown, edgeParamsUp));

    ////////////////////////////////////////////////////////////
    std::cout << "outputUp: " << outputAccessorUp[0][0] << ", " << outputAccessorUp[0][1] << ", " << outputAccessorUp[0][2] << std::endl;
    std::cout << "outputDown: " << outputAccessorDown[0][0] << ", " << outputAccessorDown[0][1] << ", " << outputAccessorDown[0][2] << std::endl;
    ////////////////////////////////////////////////////////////

    // Invoke classifier model for final output
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 6}, input);

    int insertIndex = 0;

    for (const torch::TensorAccessor<float, 2> &accessor : {outputAccessorUp, outputAccessorDown})
    {
        for (int accessorIndex = 0; accessorIndex < 3; ++accessorIndex)
        {
            input[0][insertIndex] = accessor[0][accessorIndex];
            ++insertIndex;
        }
    }

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_primaryTrackClassifierModel, {input}, output);
    torch::TensorAccessor<float, 2> outputAccessor = output.accessor<float, 2>();

    std::cout << "Track classification score: " << outputAccessor[0][0] << std::endl;
    std::cout << "------------------------" << std::endl;

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

torch::TensorAccessor<float, 2> MLPPrimaryHierarchyTool::ClassifyTrackEdge(const MLPPrimaryNetworkParams &edgeParams, 
    const MLPPrimaryNetworkParams &otherEdgeParams)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 17}, input);

    int insertIndex(this->AddToInput(0, edgeParams.GetCommonParamsForModel(), input));
    insertIndex = this->AddToInput(insertIndex, edgeParams.GetOrientationParamsForModel(), input);
    insertIndex = this->AddToInput(insertIndex, otherEdgeParams.GetOrientationParamsForModel(), input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_primaryTrackBranchModel, {input}, output);

    return output.accessor<float, 2>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPPrimaryHierarchyTool::ClassifyShower(const MLPPrimaryNetworkParams &primaryNetworkParams)
{
    std::cout << "------------------------" << std::endl;
    std::cout << "Classifying shower with " << primaryNetworkParams.m_nSpacepoints << "hits" << std::endl;

    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 9}, input);

    int insertIndex(this->AddToInput(0, primaryNetworkParams.GetCommonParamsForModel(), input));
    insertIndex = this->AddToInput(insertIndex, primaryNetworkParams.GetOrientationParamsForModel(), input);

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

    return MLPBaseHierarchyTool::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
