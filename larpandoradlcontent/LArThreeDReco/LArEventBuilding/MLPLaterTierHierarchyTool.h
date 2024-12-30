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
        float m_overshootStartDCA = -999.f;
        float m_overshootStartL = -999.f;
        float m_overshootEndDCA = -999.f;
        float m_overshootEndL = -999.f;
        float m_childConnectionDCA = -999.f;
        float m_childConnectionLRatio = -999.f;
        float m_parentConnectionPointNUpstreamHits = -999.f;
        float m_parentConnectionPointNDownstreamHits = -999.f;
        float m_parentConnectionPointNHitRatio = -999.f;
        float m_parentConnectionPointEigenvalueRatio = -999.f;
        float m_parentConnectionPointOpeningAngle = -999.f;
        float m_parentIsPOIClosestToNu = -999.f;
        float m_childIsPOIClosestToNu = -999.f;

        void Print();
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

    void GetConnectionPoint(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo, 
        const pandora::CartesianVector &parentStart, const pandora::CartesianVector &childStart, const pandora::CartesianVector &childStartDirection);

    std::pair<pandora::CartesianVector, bool> ExtrapolateChildToParent(const pandora::CartesianVector &parentPosition, const pandora::CartesianVector &childStart, 
        const pandora::CartesianVector &childStartDirection);

    bool DoesConnect(const pandora::CartesianVector &boundary1, const pandora::CartesianVector &boundary2,
        const pandora::CartesianVector &testPoint, const float buffer);

    void NormaliseNetworkParams(MLPLaterTierNetworkParams &laterTierNetworkParams);

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
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams::Print()
{
    std::cout << "ParentPOIClosestToNuVertex: " << this->m_parentIsPOIClosestToNu << std::endl;
    std::cout << "ChildPOIClosestToNuVertex: " << this->m_childIsPOIClosestToNu << std::endl;
    std::cout << "ParentTrackScore: " << this->m_parentTrackScore << std::endl; 
    std::cout << "ChildTrackScore: " << this->m_childTrackScore << std::endl;
    std::cout << "ParentNSpacepoints: " << this->m_parentNSpacepoints << std::endl;
    std::cout << "ChildNSpacepoints: " << this->m_childNSpacepoints << std::endl;
    std::cout << "Separation3D: " << this->m_separation3D << std::endl;
    std::cout << "ParentNuVertexSep: " << this->m_parentNuVertexSep << std::endl;
    std::cout << "ChildNuVertexSep: " << this->m_childNuVertexSep << std::endl;
    std::cout << "ParentEndRegionNHits: " << this->m_parentEndRegionNHits << std::endl;
    std::cout << "ParentEndRegionNParticles: " << this->m_parentEndRegionNParticles << std::endl;
    std::cout << "ParentEndRegionRToWall: " << this->m_parentEndRegionRToWall << std::endl;
    std::cout << "VertexSep: " << this->m_vertexSeparation << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_LATER_TIER_HIERARCHY_TOOL_H
