/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLPrimaryHierarchyTool.cc
 *
 *  @brief  Implementation of the DL primary hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLPrimaryHierarchyTool.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

void DLPrimaryHierarchyTool::DLPrimaryNetworkParams::AddCommonParamsToInput(int &insertIndex, LArDLHelper::TorchInput &modelInput) const
{
    modelInput[0][insertIndex] = m_nSpacepoints;
    ++insertIndex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLPrimaryHierarchyTool::DLPrimaryNetworkParams::AddOrientationParamsToInput(int &insertIndex, LArDLHelper::TorchInput &modelInput) const
{
    for (const float param : {m_nuSeparation, m_vertexRegionNHits, m_vertexRegionNParticles, m_dca, m_connectionExtrapDistance,
             m_isPOIClosestToNu, m_parentConnectionDistance, m_childConnectionDistance})
    {
        modelInput[0][insertIndex] = param;
        ++insertIndex;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DLPrimaryHierarchyTool::DLPrimaryHierarchyTool() :
    DLBaseHierarchyTool(),
    m_trainingMode(false),
    m_extrapolationStepSize(1.f),
    m_normalise(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLPrimaryHierarchyTool::Run(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pNeutrinoPfo,
    const HierarchyPfoVector &trackPfos, const HierarchyPfo &hierarchyPfo, std::vector<DLPrimaryNetworkParams> &networkParamVector, float &primaryScore)
{
    networkParamVector.clear();
    primaryScore = std::numeric_limits<float>::lowest();

    this->SetDetectorBoundaries();

    const bool isTrack(hierarchyPfo.GetPfo()->GetParticleId() == 13);
    std::vector<bool> orientationVector(isTrack ? std::vector<bool>({true, false}) : std::vector<bool>({true}));

    for (const bool useUpstream : orientationVector)
    {
        // Set network params
        DLPrimaryNetworkParams primaryNetworkParams;

        const StatusCode statusCode(this->CalculateNetworkVariables(pAlgorithm, hierarchyPfo, pNeutrinoPfo, trackPfos, useUpstream, primaryNetworkParams));

        if (statusCode != STATUS_CODE_SUCCESS)
            return statusCode;

        networkParamVector.emplace_back(primaryNetworkParams);
    }

    // Now run the model!
    if (!m_trainingMode)
    {
        if (isTrack)
            primaryScore = this->ClassifyTrack(networkParamVector.at(0), networkParamVector.at(1));
        else
            primaryScore = this->ClassifyShower(networkParamVector.at(0));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLPrimaryHierarchyTool::CalculateNetworkVariables(const Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo,
    const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos, const bool useUpstream,
    DLPrimaryNetworkParams &primaryNetworkParams) const
{
    // Pick out neutrino vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;

    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(
        pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    // Pick out the correct extremal point
    const ExtremalPoint &particlePoint(useUpstream ? hierarchyPfo.GetUpstreamPoint() : hierarchyPfo.GetDownstreamPoint());

    // Check that we could actually calculate pfo vertex/direction, return if not
    if (!particlePoint.IsSet())
    {
        return STATUS_CODE_NOT_INITIALIZED;
    }

    // Set primaryNetworkParams
    primaryNetworkParams.m_isPOIClosestToNu = useUpstream ? 1.f : 0.f;
    primaryNetworkParams.m_nSpacepoints = LArPfoHelper::GetNumberOfThreeDHits(hierarchyPfo.GetPfo());
    primaryNetworkParams.m_nuSeparation = (particlePoint.GetPosition() - nuVertex).GetMagnitude();
    this->SetVertexRegionParams(pAlgorithm, hierarchyPfo.GetPfo(), particlePoint.GetPosition(), primaryNetworkParams);
    this->SetConnectionParams(particlePoint, nuVertex, primaryNetworkParams);
    this->SetContextParams(hierarchyPfo.GetPfo(), particlePoint, nuVertex, trackPfos, primaryNetworkParams);

    // Normalise
    if (m_normalise)
        this->NormaliseNetworkParams(primaryNetworkParams);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLPrimaryHierarchyTool::SetVertexRegionParams(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    const CartesianVector &particleVertex, DLPrimaryNetworkParams &primaryNetworkParams) const
{
    std::pair<float, float> vertexRegionParams(this->GetParticleInfoAboutPfoPosition(pAlgorithm, pPfo, particleVertex));

    primaryNetworkParams.m_vertexRegionNHits = vertexRegionParams.first;
    primaryNetworkParams.m_vertexRegionNParticles = vertexRegionParams.second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLPrimaryHierarchyTool::SetConnectionParams(
    const ExtremalPoint &particlePoint, const CartesianVector &nuVertex, DLPrimaryNetworkParams &primaryNetworkParams) const
{
    // Extrapolate particle to nuVertex
    const float extrapDistance((nuVertex - particlePoint.GetPosition()).GetDotProduct(particlePoint.GetDirection()));
    const CartesianVector extrapolationPoint(particlePoint.GetPosition() + (particlePoint.GetDirection() * extrapDistance));

    primaryNetworkParams.m_dca = (nuVertex - extrapolationPoint).GetMagnitude();
    primaryNetworkParams.m_connectionExtrapDistance = (extrapDistance * (-1.f)); // backwards extrap. should be +ve
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLPrimaryHierarchyTool::SetContextParams(const ParticleFlowObject *const pPfo, const ExtremalPoint &particlePoint,
    const CartesianVector &nuVertex, const HierarchyPfoVector &trackPfos, DLPrimaryNetworkParams &primaryNetworkParams) const
{
    // What is our nuVertex separation
    const float nuVertexSepSq((particlePoint.GetPosition() - nuVertex).GetMagnitudeSquared());

    bool found(false);
    float parentConnectionDistance(std::numeric_limits<float>::max());
    float childConnectionDistance(std::numeric_limits<float>::max());

    // Look at potential parents
    for (const HierarchyPfo &hierarchyTrackPfo : trackPfos)
    {
        if (hierarchyTrackPfo == pPfo)
            continue;

        // Upstream position should be closer to nu vertex
        const float thisNuVertexSepSq((hierarchyTrackPfo.GetDownstreamPoint().GetPosition() - nuVertex).GetMagnitudeSquared());

        if (thisNuVertexSepSq > nuVertexSepSq)
            continue;

        float thisParentConnectionDistance(std::numeric_limits<float>::lowest());
        float thisChildConnectionDistance(std::numeric_limits<float>::lowest());

        this->CalculateConnectionDistances(hierarchyTrackPfo.GetDownstreamPoint(), particlePoint, thisParentConnectionDistance, thisChildConnectionDistance);

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

void DLPrimaryHierarchyTool::CalculateConnectionDistances(
    const ExtremalPoint &parentPoint, const ExtremalPoint &childPoint, float &parentConnectionDistance, float &childConnectionDistance) const
{
    // Loop things
    float smallestT(std::numeric_limits<float>::max());
    bool isGettingCloser(true);
    bool found(false);
    CartesianVector connectionPoint(std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest());

    // start the seed
    CartesianVector extrapolatedPoint(childPoint.GetPosition());

    while (isGettingCloser)
    {
        isGettingCloser = false;
        extrapolatedPoint = extrapolatedPoint + (childPoint.GetDirection() * ((-1.f) * m_extrapolationStepSize));

        if (!LArGeometryHelper::IsInDetector(m_detectorBoundaries, extrapolatedPoint))
            break;

        const float parentT(
            (parentPoint.GetDirection() * (-1.f)).GetCrossProduct((extrapolatedPoint - parentPoint.GetPosition())).GetMagnitudeSquared());

        if (parentT < smallestT)
        {
            smallestT = parentT;
            connectionPoint = extrapolatedPoint;
            isGettingCloser = true;
            found = true;
        }
    }

    // Distance from parent end to connection point - need to turn parent direction around
    parentConnectionDistance = found ? (parentPoint.GetDirection() * (-1.f)).GetDotProduct(connectionPoint - parentPoint.GetPosition())
                                     : std::numeric_limits<float>::lowest();
    // Distance from child start to connection point
    childConnectionDistance = found ? (childPoint.GetPosition() - connectionPoint).GetMagnitude() : std::numeric_limits<float>::lowest();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLPrimaryHierarchyTool::NormaliseNetworkParams(DLPrimaryNetworkParams &primaryNetworkParams) const
{
    this->NormaliseNetworkParam(m_normLimits.m_nSpacepointsMin, m_normLimits.m_nSpacepointsMax, primaryNetworkParams.m_nSpacepoints);
    this->NormaliseNetworkParam(m_normLimits.m_nuSeparationMin, m_normLimits.m_nuSeparationMax, primaryNetworkParams.m_nuSeparation);
    this->NormaliseNetworkParam(m_normLimits.m_vertexRegionNHitsMin, m_normLimits.m_vertexRegionNHitsMax, primaryNetworkParams.m_vertexRegionNHits);
    this->NormaliseNetworkParam(
        m_normLimits.m_vertexRegionNParticlesMin, m_normLimits.m_vertexRegionNParticlesMax, primaryNetworkParams.m_vertexRegionNParticles);
    this->NormaliseNetworkParam(m_normLimits.m_dcaMin, m_normLimits.m_dcaMax, primaryNetworkParams.m_dca);
    this->NormaliseNetworkParam(m_normLimits.m_connectionExtrapDistanceMin, m_normLimits.m_connectionExtrapDistanceMax,
        primaryNetworkParams.m_connectionExtrapDistance);
    this->NormaliseNetworkParam(m_normLimits.m_parentConnectionDistanceMin, m_normLimits.m_parentConnectionDistanceMax,
        primaryNetworkParams.m_parentConnectionDistance);
    this->NormaliseNetworkParam(m_normLimits.m_childConnectionDistanceMin, m_normLimits.m_childConnectionDistanceMax,
        primaryNetworkParams.m_childConnectionDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLPrimaryHierarchyTool::ClassifyTrack(const DLPrimaryNetworkParams &edgeParamsUp, const DLPrimaryNetworkParams &edgeParamsDown)
{
    // Invoke branch model for each edge
    const FloatVector outputUp(this->ClassifyTrackEdge(edgeParamsUp, edgeParamsDown));
    const FloatVector outputDown(this->ClassifyTrackEdge(edgeParamsDown, edgeParamsUp));

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
    LArDLHelper::Forward(m_primaryTrackClassifierModel, {input}, output);
    torch::TensorAccessor<float, 2> outputAccessor = output.accessor<float, 2>();

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector DLPrimaryHierarchyTool::ClassifyTrackEdge(const DLPrimaryNetworkParams &edgeParams, const DLPrimaryNetworkParams &otherEdgeParams)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 17}, input);

    // Fill our torch vector
    int insertIndex(0);
    edgeParams.AddCommonParamsToInput(insertIndex, input);
    edgeParams.AddOrientationParamsToInput(insertIndex, input);
    otherEdgeParams.AddOrientationParamsToInput(insertIndex, input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_primaryTrackBranchModel, {input}, output);

    torch::TensorAccessor<float, 2> outputAccessor(output.accessor<float, 2>());

    return {outputAccessor[0][0], outputAccessor[0][1], outputAccessor[0][2]};
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLPrimaryHierarchyTool::ClassifyShower(const DLPrimaryNetworkParams &primaryNetworkParams)
{
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 9}, input);

    // Fill our torch vector
    int insertIndex(0);
    primaryNetworkParams.AddCommonParamsToInput(insertIndex, input);
    primaryNetworkParams.AddOrientationParamsToInput(insertIndex, input);

    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_primaryShowerClassifierModel, {input}, output);
    auto outputAccessor = output.accessor<float, 2>();

    return outputAccessor[0][0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLPrimaryHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ExtrapolationStepSize", m_extrapolationStepSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Normalise", m_normalise));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMin", m_normLimits.m_nSpacepointsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMax", m_normLimits.m_nSpacepointsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuSeparationMin", m_normLimits.m_nuSeparationMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuSeparationMax", m_normLimits.m_nuSeparationMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "VertexRegionNHitsMin", m_normLimits.m_vertexRegionNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "VertexRegionNHitsMax", m_normLimits.m_vertexRegionNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "VertexRegionNParticlesMin", m_normLimits.m_vertexRegionNParticlesMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "VertexRegionNParticlesMax", m_normLimits.m_vertexRegionNParticlesMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DCAMin", m_normLimits.m_dcaMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DCAMax", m_normLimits.m_dcaMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ConnectionExtrapDistanceMin", m_normLimits.m_connectionExtrapDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ConnectionExtrapDistanceMax", m_normLimits.m_connectionExtrapDistanceMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentConnectionDistanceMin", m_normLimits.m_parentConnectionDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ParentConnectionDistanceMax", m_normLimits.m_parentConnectionDistanceMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ChildConnectionDistanceMin", m_normLimits.m_childConnectionDistanceMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ChildConnectionDistanceMax", m_normLimits.m_childConnectionDistanceMax));

    if (!m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryTrackBranchModelName", m_primaryTrackBranchModelName));
        m_primaryTrackBranchModelName = LArFileHelper::FindFileInPath(m_primaryTrackBranchModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_primaryTrackBranchModelName, m_primaryTrackBranchModel));

        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryTrackClassifierModelName", m_primaryTrackClassifierModelName));
        m_primaryTrackClassifierModelName = LArFileHelper::FindFileInPath(m_primaryTrackClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_primaryTrackClassifierModelName, m_primaryTrackClassifierModel));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            XmlHelper::ReadValue(xmlHandle, "PrimaryShowerClassifierModelName", m_primaryShowerClassifierModelName));
        m_primaryShowerClassifierModelName = LArFileHelper::FindFileInPath(m_primaryShowerClassifierModelName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_primaryShowerClassifierModelName, m_primaryShowerClassifierModel));
    }

    return DLBaseHierarchyTool::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
