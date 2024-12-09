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
        float m_closestParentL = -999.f;
        float m_closestParentT = -999.f;
    };

    /**
     *  @brief  Default constructor
     */
    MLPPrimaryHierarchyTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo, const pandora::ParticleFlowObject *const pNeutrinoPfo, 
        const bool useUpstream);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void SetNSpacepoints(const pandora::ParticleFlowObject *const pPfo, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void SetNuVertexSep(const pandora::CartesianVector &particleVertex, const pandora::ParticleFlowObject *const pNeutrinoPfo,
        MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void SetVertexRegionParams(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo, 
        const pandora::CartesianVector &particleVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void NormaliseNetworkParams(MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void NormaliseNetworkParam(const float minLimit, const float maxLimit, float &primaryNetworkParam) const;

    pandora::StringVector m_pfoListNames;
    float m_nSpacepointsMin;
    float m_nSpacepointsMax;
    float m_nuSeparationMin;
    float m_nuSeparationMax;
    float m_vertexRegionRadius;
    float m_vertexRegionNHitsMin;
    float m_vertexRegionNHitsMax;
    float m_vertexRegionNParticlesMin;
    float m_vertexRegionNParticlesMax;
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_PRIMARY_HIERARCHY_TOOL_H
