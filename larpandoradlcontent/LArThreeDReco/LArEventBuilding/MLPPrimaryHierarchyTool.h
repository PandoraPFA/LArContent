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

namespace lar_dl_content
{

/**
 *   @brief  MLPPrimaryHierarchyTool to calculate variables related to the initial shower region
 */
class MLPPrimaryHierarchyTool : public pandora::AlgorithmTool
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
    };

    /**
     *  @brief  Default constructor
     */
    MLPPrimaryHierarchyTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pNeutrinoPfo, 
        const HierarchyPfoMap &trackPfos, HierarchyPfo &hierarchyPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void SetDetectorBoundaries();

    pandora::StatusCode CalculateNetworkVariables(const pandora::Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo, 
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const bool useUpstream, 
        MLPPrimaryNetworkParams &primaryNetworkParams);

    void SetNSpacepoints(const pandora::ParticleFlowObject *const pPfo, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void SetNuVertexSep(const pandora::CartesianVector &particleVertex, const pandora::CartesianVector &nuVertex,
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

    bool IsInFV(const pandora::CartesianVector &position) const;

    void NormaliseNetworkParams(MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void NormaliseNetworkParam(const float minLimit, const float maxLimit, float &primaryNetworkParam) const;

    float ClassifyTrack(const MLPPrimaryNetworkParams &primaryNetworkParamsUp, const MLPPrimaryNetworkParams &primaryNetworkParamsDown);

    float ClassifyShower(const MLPPrimaryNetworkParams &primaryNetworkParams);

    std::string m_primaryTrackBranchModelName;
    std::string m_primaryTrackClassifierModelName;
    std::string m_primaryShowerClassifierModelName;

    LArDLHelper::TorchModel m_primaryTrackBranchModel;
    LArDLHelper::TorchModel m_primaryTrackClassifierModel;
    LArDLHelper::TorchModel m_primaryShowerClassifierModel;

    float m_detectorMinX;
    float m_detectorMaxX;
    float m_detectorMinY;
    float m_detectorMaxY;
    float m_detectorMinZ;
    float m_detectorMaxZ;
    pandora::StringVector m_pfoListNames;
    float m_vertexRegionRadius;
    float m_extrapolationStepSize;
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

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_PRIMARY_HIERARCHY_TOOL_H
