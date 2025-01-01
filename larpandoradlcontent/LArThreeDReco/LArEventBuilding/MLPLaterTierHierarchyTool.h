/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.h
 *
 *  @brief  Header file for the MLP later hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_MLP_LATER_TIER_HIERARCHY_TOOL_H
#define LAR_MLP_LATER_TIER_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPBaseHierarchyTool.h"

namespace lar_dl_content
{

/**
 *   @brief  MLPLaterTierHierarchyTool to calculate variables related to the initial shower region
 */
class MLPLaterTierHierarchyTool : public MLPBaseHierarchyTool
{
public:
    struct MLPLaterTierNetworkParams  // all floats because they'll be normalised
    {
        float m_parentTrackScore = -999.f;
        float m_childTrackScore = -999.f;
        float m_parentNSpacepoints = -999.f;
        float m_childNSpacepoints = -999.f;
        float m_separation3D = -999.f;
        float m_parentNuVertexSep = -999.f;
        float m_childNuVertexSep = -999.f;
        float m_parentEndRegionNHits = -999.f;
        float m_parentEndRegionNParticles = -999.f;
        float m_parentEndRegionRToWall = -999.f;
        float m_vertexSeparation = -999.f;
        float m_doesChildConnect = -999.f;
        pandora::CartesianVector m_connectionPoint = pandora::CartesianVector(-999.f, -999.f, -999.f);
        pandora::CartesianVector m_connectionDirection = pandora::CartesianVector(-999.f, -999.f, -999.f);
        float m_overshootStartDCA = -999.f;
        float m_overshootStartL = -999.f;
        float m_overshootEndDCA = -999.f;
        float m_overshootEndL = -999.f;
        float m_childCPDCA = -999.f;
        float m_childCPExtrapDistance = -999.f;
        float m_childCPLRatio = -999.f;
        float m_parentCPNUpstreamHits = -999.f;
        float m_parentCPNDownstreamHits = -999.f;
        float m_parentCPNHitRatio = -999.f;
        float m_parentCPEigenvalueRatio = -999.f;
        float m_parentCPOpeningAngle = -999.f;
        float m_parentIsPOIClosestToNu = -999.f;
        float m_childIsPOIClosestToNu = -999.f;

        void Print();

        pandora::FloatVector GetCommonParamsForModel() const;

        pandora::FloatVector GetOrientationParamsForModel() const;
    };

    /**
     *  @brief  Default constructor
     */
    MLPLaterTierHierarchyTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pNeutrinoPfo, 
        HierarchyPfo &parentHierarchyPfo, HierarchyPfo &childHierarchyPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool IsShowerVertexUpstream(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo);

    pandora::StatusCode CalculateNetworkVariables(const pandora::Algorithm *const pAlgorithm, 
        const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const bool useUpstreamForParent, const bool useUpstreamForChild, 
        MLPLaterTierNetworkParams &laterTierNetworkParams);

    std::pair<float, float> GetTrackScoreParams(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo);
    std::pair<float, float> GetNSpacepointsParams(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo);
    float GetSeparation3D(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo);

    void SetCommonParams(const std::pair<float, float> &trackScoreParams, const std::pair<float, float> &nSpacepointsParams, 
        const float separation3D, MLPLaterTierNetworkParams &laterTierNetworkParams);

    void SetVertexParams(const pandora::CartesianVector &nuVertex, const pandora::CartesianVector &parentStart, 
        const pandora::CartesianVector &parentEnd, const pandora::CartesianVector &childStart, 
        MLPLaterTierNetworkParams &laterTierNetworkParams);

    void SetEndRegionParams(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pParentPfo, 
        const pandora::CartesianVector &parentEnd, MLPLaterTierNetworkParams &laterTierNetworkParams) const;

    void SetEndRegionRToWall(const pandora::CartesianVector &parentEnd, MLPLaterTierNetworkParams &laterTierNetworkParams) const;

    void SetConnectionParams(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo, 
        const pandora::CartesianVector &parentStart, const pandora::CartesianVector &childStart, 
        const pandora::CartesianVector &childStartDirection, MLPLaterTierNetworkParams &laterTierNetworkParams);

    std::pair<pandora::CartesianVector, bool> ExtrapolateChildToParent(const pandora::CartesianVector &parentPosition, 
        const pandora::CartesianVector &childStart, const pandora::CartesianVector &childStartDirection);

    bool DoesConnect(const pandora::CartesianVector &boundary1, const pandora::CartesianVector &boundary2,
        const pandora::CartesianVector &testPoint, const float buffer);

    void SetOvershootParams(const pandora::CartesianVector &parentStart, const pandora::CartesianVector &parentStartDirection, 
        const pandora::CartesianVector &parentEnd, const pandora::CartesianVector &parentEndDirection, const pandora::CartesianVector &childStart, 
        const pandora::CartesianVector &childStartDirection, MLPLaterTierNetworkParams &laterTierNetworkParams);

    void SetParentConnectionPointVars(const HierarchyPfo &parentHierarchyPfo, MLPLaterTierNetworkParams &laterTierNetworkParams);

    void NormaliseNetworkParams(MLPLaterTierNetworkParams &laterTierNetworkParams);

    float ClassifyTrackTrack(const MLPLaterTierNetworkParams &edgeParamsUpUp, const MLPLaterTierNetworkParams &edgeParamsUpDown, 
        const MLPLaterTierNetworkParams &edgeParamsDownUp, const MLPLaterTierNetworkParams &edgeParamsDownDown);

    pandora::FloatVector ClassifyTrackTrackEdge(const MLPLaterTierNetworkParams &edgeParams, 
        const MLPLaterTierNetworkParams &otherEdgeParams1, const MLPLaterTierNetworkParams &otherEdgeParams2, 
        const MLPLaterTierNetworkParams &otherEdgeParams3);

    float ClassifyTrackShower(const MLPLaterTierNetworkParams &edgeParamsUp, const MLPLaterTierNetworkParams &edgeParamsDown);

    pandora::FloatVector ClassifyTrackShowerEdge(const MLPLaterTierNetworkParams &edgeParams, 
        const MLPLaterTierNetworkParams &otherEdgeParams);

    // For model
    std::string m_trackTrackBranchModelName;
    std::string m_trackTrackClassifierModelName;
    std::string m_trackShowerBranchModelName;
    std::string m_trackShowerClassifierModelName;
    LArDLHelper::TorchModel m_trackTrackBranchModel;
    LArDLHelper::TorchModel m_trackTrackClassifierModel;
    LArDLHelper::TorchModel m_trackShowerBranchModel;
    LArDLHelper::TorchModel m_trackShowerClassifierModel;
    // For tool
    float m_trajectoryStepSize;
    float m_connectionBuffer;
    float m_searchRegion;
    // For normalisation
    float m_trackScoreMin;
    float m_trackScoreMax;
    float m_nSpacepointsMin;
    float m_nSpacepointsMax;
    float m_separation3DMin;
    float m_separation3DMax;
    float m_nuVertexSepMin;
    float m_nuVertexSepMax;
    float m_parentEndRegionNHitsMin;
    float m_parentEndRegionNHitsMax;
    float m_parentEndRegionNParticlesMin;
    float m_parentEndRegionNParticlesMax;
    float m_parentEndRegionRToWallMin;
    float m_parentEndRegionRToWallMax;
    float m_vertexSepMin;
    float m_vertexSepMax;
    float m_doesChildConnectMin;
    float m_doesChildConnectMax;
    float m_overshootDCAMin;
    float m_overshootDCAMax;
    float m_overshootLMin;
    float m_overshootLMax;
    float m_childCPDCAMin;
    float m_childCPDCAMax;
    float m_childCPExtrapDistanceMin;
    float m_childCPExtrapDistanceMax;
    float m_childCPLRatioMin;
    float m_childCPLRatioMax;
    float m_parentCPNHitsMin;
    float m_parentCPNHitsMax;
    float m_parentCPNHitRatioMin;
    float m_parentCPNHitRatioMax;
    float m_parentCPEigenvalueRatioMin;
    float m_parentCPEigenvalueRatioMax;
    float m_parentCPOpeningAngleMin;
    float m_parentCPOpeningAngleMax;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams::Print()
{
    std::cout << "ParentPOIClosestToNuVertex: " << m_parentIsPOIClosestToNu << std::endl;
    std::cout << "ChildPOIClosestToNuVertex: " << m_childIsPOIClosestToNu << std::endl;
    std::cout << "ParentTrackScore: " << m_parentTrackScore << std::endl; 
    std::cout << "ChildTrackScore: " << m_childTrackScore << std::endl;
    std::cout << "ParentNSpacepoints: " << m_parentNSpacepoints << std::endl;
    std::cout << "ChildNSpacepoints: " << m_childNSpacepoints << std::endl;
    std::cout << "Separation3D: " << m_separation3D << std::endl;
    std::cout << "ParentNuVertexSep: " << m_parentNuVertexSep << std::endl;
    std::cout << "ChildNuVertexSep: " << m_childNuVertexSep << std::endl;
    std::cout << "ParentEndRegionNHits: " << m_parentEndRegionNHits << std::endl;
    std::cout << "ParentEndRegionNParticles: " << m_parentEndRegionNParticles << std::endl;
    std::cout << "ParentEndRegionRToWall: " << m_parentEndRegionRToWall << std::endl;
    std::cout << "VertexSep: " << m_vertexSeparation << std::endl;
    std::cout << "DoesChildConnect: " << m_doesChildConnect << std::endl;
    std::cout << "OvershootStartDCA: " << m_overshootStartDCA << std::endl;
    std::cout << "OvershootStartL: " << m_overshootStartL << std::endl;
    std::cout << "OvershootEndDCA: " << m_overshootEndDCA << std::endl;
    std::cout << "OvershootEndL: " << m_overshootEndL << std::endl;
    std::cout << "ChildCPDCA: " << m_childCPDCA << std::endl;
    std::cout << "ChildCPExtrapDistance: " << m_childCPExtrapDistance << std::endl;
    std::cout << "ChildCPLRatio: " << m_childCPLRatio << std::endl;
    std::cout << "ParentCPNUpstreamHits: " << m_parentCPNUpstreamHits << std::endl;
    std::cout << "ParentCPNDownstreamHits: " << m_parentCPNDownstreamHits << std::endl;
    std::cout << "ParentCPNHitRatio: " << m_parentCPNHitRatio << std::endl;
    std::cout << "ParentCPEigenvalueRatio: " << m_parentCPEigenvalueRatio << std::endl;
    std::cout << "ParentCPOpeningAngle: " << m_parentCPOpeningAngle << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams::GetCommonParamsForModel() const
{
    return {m_parentTrackScore, m_childTrackScore, m_parentNSpacepoints, m_childNSpacepoints, m_separation3D};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams::GetOrientationParamsForModel() const
{
    return {m_parentNuVertexSep, m_childNuVertexSep, m_parentEndRegionNHits, m_parentEndRegionNParticles, m_parentEndRegionRToWall, 
            m_vertexSeparation, m_doesChildConnect, m_overshootStartDCA, m_overshootStartL, m_overshootEndDCA, m_overshootEndL, 
            m_childCPDCA, m_childCPExtrapDistance, m_childCPLRatio, m_parentCPNUpstreamHits, m_parentCPNDownstreamHits, 
            m_parentCPNHitRatio, m_parentCPEigenvalueRatio, m_parentCPOpeningAngle, m_parentIsPOIClosestToNu, m_childIsPOIClosestToNu};
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_LATER_TIER_HIERARCHY_TOOL_H
