/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h
 *
 *  @brief  Header file for the shower start refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_REFINEMENT_ALGORITHM_H
#define LAR_SHOWER_START_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

#include "TMVA/Reader.h"

namespace lar_content
{

class ProtoShower;
class ShowerStartRefinementBaseTool;

class ShowerStartRefinementAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    ShowerStartRefinementAlgorithm();
    ~ShowerStartRefinementAlgorithm();

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
    float m_electronFraction;

    pandora::CaloHitList GetAllHitsOfType(const pandora::HitType hitType);
    pandora::CaloHitList GetXIntervalHitsOfType(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType);

    void FillOwnershipMaps();

    bool IsElectronPathway(const pandora::CaloHitList &hitsToAdd);
    void AddElectronPathway(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CaloHitList &pathwayHitList);

    void SetElectronTreeMetadata(const pandora::ParticleFlowObject *const pShowerPfo, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);
    void SetElectronMetadata(const pandora::CartesianVector &nuVertexPosition, const pandora::ParticleFlowObject *const pShowerPfo);
    void SetGammaVertex(const pandora::CartesianVector &showerVertex, const pandora::ParticleFlowObject *const pShowerPfo);
    bool IsTrack(const ProtoShower &protoShower);
    bool TMVAIsElectron(LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables, const pandora::ParticleFlowObject *const pShowerPfo, const bool alterMetadata);
    bool TMVAIsGamma(LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables, const pandora::ParticleFlowObject *const pShowerPfo);
    void RemoveConnectionPathway(const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower);

    void FillGammaHitMap();
    void FillElectronHitMap();
    bool IsElectron(const pandora::ParticleFlowObject *const pPfo) const;
    bool IsGamma(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &nuVertexPosition) const;
    void FillTree(const std::string &treeName, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables);

    bool m_createTrainingTrees;
    bool m_hybridMode;
    float m_electronTMVACut;
    float m_gammaTMVACut;
    LArConnectionPathwayHelper::ElectronTreeVariables m_electronTreeVariables;

private:
    pandora::StatusCode Run();
    void InitialiseElectronTrees();
    void FillPfoVector(pandora::PfoVector &pfoVector);
    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &neutrinoVertex);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;
    std::string m_neutrinoVertexListName;

    typedef std::vector<ShowerStartRefinementBaseTool *> ShowerStartRefinementToolVector;
    ShowerStartRefinementToolVector m_algorithmToolVector;

    float m_minElectronCompleteness;
    float m_minElectronPurity;
    float m_minGammaCompleteness;
    float m_thresholdSignalGammaDisplacement;

    TMVA::Reader m_TMVAReader;
    LArConnectionPathwayHelper::ElectronTreeVariables m_TMVAElectronTreeVariables;

    std::map<const pandora::MCParticle*, pandora::CaloHitList> m_gammaHitMap;
    std::map<const pandora::MCParticle*, pandora::CaloHitList> m_electronHitMap;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
