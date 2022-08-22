/**
 *  @file   larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementAlgorithm.h
 *
 *  @brief  Header file for the shower start refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HYBRID_SHOWER_START_REFINEMENT_ALGORITHM_H
#define LAR_HYBRID_SHOWER_START_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

namespace lar_content
{

class ProtoShower;
class HybridShowerStartRefinementBaseTool;

class HybridShowerStartRefinementAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    HybridShowerStartRefinementAlgorithm();

    typedef std::unordered_map<const pandora::MCParticle*, pandora::CaloHitList> MCParticleToHitListMap;

    typedef std::map<int, int> DeviationAngleMap;
    DeviationAngleMap m_thetaMapU;
    DeviationAngleMap m_thetaMapV;
    DeviationAngleMap m_thetaMapW;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::unordered_map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    HitToClusterMap m_hitToClusterMapU;
    HitToClusterMap m_hitToClusterMapV;
    HitToClusterMap m_hitToClusterMapW;
    ClusterToPfoMap m_clusterToPfoMapU;
    ClusterToPfoMap m_clusterToPfoMapV;
    ClusterToPfoMap m_clusterToPfoMapW;

    pandora::PfoList m_deletedPfos;
    float m_binSize;

    pandora::CaloHitList GetAllHitsOfType(const pandora::HitType hitType);
    pandora::CaloHitList GetXIntervalHitsOfType(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType);

    void FillOwnershipMaps();

    void AddElectronPathway(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CaloHitList &pathwayHitList);

    void SetElectronMetadata(const pandora::CartesianVector &nuVertexPosition, const pandora::ParticleFlowObject *const pShowerPfo);
    //void RemoveConnectionPathway(const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower);
    void RemoveConnectionPathway(const pandora::ParticleFlowObject *const pGammaPfo, const ElectronProtoShowerVector &protoShowerVectorU, 
        const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW);

    bool IsShower(const pandora::MCParticle *const pMCParticle);
    const pandora::Cluster *CreateCluster(const pandora::MCParticle *const pMCParticle, const pandora::CaloHitList &caloHitList, const pandora::HitType hitType);

    void FillGammaHitMap();
    void FillElectronHitMap();
    bool IsElectron(const pandora::ParticleFlowObject *const pPfo) const;
    bool IsGamma(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &nuVertexPosition) const;

private:
    pandora::StatusCode Run();
    void FillPfoVector(pandora::PfoVector &pfoVector);
    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &neutrinoVertex);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;
    std::string m_neutrinoVertexListName;

    typedef std::vector<HybridShowerStartRefinementBaseTool *> ShowerStartRefinementToolVector;
    ShowerStartRefinementToolVector m_algorithmToolVector;

    float m_minElectronCompleteness;
    float m_minElectronPurity;
    float m_minGammaCompleteness;
    float m_thresholdSignalGammaDisplacement;

    std::map<const pandora::MCParticle*, pandora::CaloHitList> m_gammaHitMap;
    std::map<const pandora::MCParticle*, pandora::CaloHitList> m_electronHitMap;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
