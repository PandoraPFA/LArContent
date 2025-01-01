/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h
 *
 *  @brief  Header file for the MLP primary hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_MLP_PRIMARY_HIERARCHY_TOOL_H
#define LAR_MLP_PRIMARY_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPBaseHierarchyTool.h"

namespace lar_dl_content
{

/**
 *   @brief  MLPPrimaryHierarchyTool to calculate variables related to the initial shower region
 */
class MLPPrimaryHierarchyTool : public MLPBaseHierarchyTool
{
public:
    struct MLPPrimaryNetworkParams  // all floats because they'll be normalised
    {
        float m_nSpacepoints = -999.f;
        float m_nuSeparation = -999.f;
        float m_vertexRegionNHits = -999.f;
        float m_vertexRegionNParticles = -999.f;
        float m_dca = -999.f;
        float m_connectionExtrapDistance = -999.f;
        float m_isPOIClosestToNu = -999.f;
        float m_parentConnectionDistance = -999.f;
        float m_childConnectionDistance = -999.f;

        void Print() const;

        pandora::FloatVector GetCommonParamsForModel() const;

        pandora::FloatVector GetOrientationParamsForModel() const;
    };

    /**
     *  @brief  Default constructor
     */
    MLPPrimaryHierarchyTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pNeutrinoPfo, 
        const HierarchyPfoMap &trackPfos, const HierarchyPfo &hierarchyPfo, float &primaryScore);

private:

    pandora::StatusCode CalculateNetworkVariables(const pandora::Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo, 
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const bool useUpstream, 
        MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void SetVertexRegionParams(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo, 
        const pandora::CartesianVector &particleVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void SetConnectionParams(const pandora::CartesianVector &particleVertex, const pandora::CartesianVector &particleDirection, 
        const pandora::CartesianVector &nuVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void SetContextParams(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &particleVertex, 
        const pandora::CartesianVector &particleDirection, const pandora::CartesianVector &nuVertex, 
        const HierarchyPfoMap &trackPfos, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void CalculateConnectionDistances(const pandora::CartesianVector &parentVertex, const pandora::CartesianVector &parentDirection, 
        const pandora::CartesianVector &childVertex, const pandora::CartesianVector &childDirection, 
        float &trainingCutL, float &trainingCutT) const;

    void NormaliseNetworkParams(MLPPrimaryNetworkParams &primaryNetworkParams) const;

    float ClassifyTrack(const MLPPrimaryNetworkParams &edgeParamsUp, const MLPPrimaryNetworkParams &edgeParamsDown);

    pandora::FloatVector ClassifyTrackEdge(const MLPPrimaryNetworkParams &edgeParams, 
        const MLPPrimaryNetworkParams &otherEdgeParams);

    float ClassifyShower(const MLPPrimaryNetworkParams &primaryNetworkParams);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // For model
    std::string m_primaryTrackBranchModelName;
    std::string m_primaryTrackClassifierModelName;
    std::string m_primaryShowerClassifierModelName;
    LArDLHelper::TorchModel m_primaryTrackBranchModel;
    LArDLHelper::TorchModel m_primaryTrackClassifierModel;
    LArDLHelper::TorchModel m_primaryShowerClassifierModel;
    // For tool
    float m_extrapolationStepSize;
    // For normalisation
    float m_nSpacepointsMin;
    float m_nSpacepointsMax;
    float m_nuSeparationMin;
    float m_nuSeparationMax;
    float m_vertexRegionNHitsMin;
    float m_vertexRegionNHitsMax;
    float m_vertexRegionNParticlesMin;
    float m_vertexRegionNParticlesMax;
    float m_dcaMin;
    float m_dcaMax;
    float m_connectionExtrapDistanceMin;
    float m_connectionExtrapDistanceMax;
    float m_parentConnectionDistanceMin;
    float m_parentConnectionDistanceMax;
    float m_childConnectionDistanceMin;
    float m_childConnectionDistanceMax;
};

//------------------------------------------------------------------------------------------------------------------------------------------

 inline void MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams::Print() const
{
    std::cout << "IsPOIClosestToNu: " << m_isPOIClosestToNu << std::endl;
    std::cout << "NSpacepoints: " << m_nSpacepoints << std::endl;
    std::cout << "NuVertexSep: " << m_nuSeparation << std::endl;
    std::cout << "StartRegionNHits: " << m_vertexRegionNHits << std::endl;
    std::cout << "StartRegionNParticles: " << m_vertexRegionNParticles << std::endl;
    std::cout << "DCA: " << m_dca << std::endl;
    std::cout << "ExtrapDistance: " << m_connectionExtrapDistance << std::endl;
    std::cout << "parentConnectionDistance: " << m_parentConnectionDistance << std::endl;
    std::cout << "childConnectionDistance: " << m_childConnectionDistance << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams::GetCommonParamsForModel() const
{
    return {m_nSpacepoints};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams::GetOrientationParamsForModel() const
{
    return {m_nuSeparation, m_vertexRegionNHits, m_vertexRegionNParticles, m_dca, m_connectionExtrapDistance, 
            m_isPOIClosestToNu, m_parentConnectionDistance, m_childConnectionDistance};
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_PRIMARY_HIERARCHY_TOOL_H
