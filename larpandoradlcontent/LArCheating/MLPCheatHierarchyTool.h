/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPCheatHierarchyTool.h
 *
 *  @brief  Header file for the cheat hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_MLP_CHEAT_HIERARCHY_TOOL_H
#define LAR_MLP_CHEAT_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

using namespace lar_content;

namespace lar_dl_content
{

/**
 *   @brief  MLPCheatHierarchyTool to calculate variables related to the initial shower region
 */
class MLPCheatHierarchyTool : public pandora::AlgorithmTool
{
public:

    typedef std::map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> PfoToMCParticleMap;
    typedef std::map<const pandora::ParticleFlowObject*, const pandora::ParticleFlowObject*> PfoToPfoMap;
    
    /**
     *  @brief  Default constructor
     */
    MLPCheatHierarchyTool();

    pandora::StatusCode Run(const PfoToMCParticleMap &pfoToMCParticleMap, const PfoToPfoMap &childToParentPfoMap,
        const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo, bool &isTrueLink, bool &trueParentOrientation,
        bool &trueChildOrientation);

    pandora::StatusCode Run(const PfoToMCParticleMap &pfoToMCParticleMap, const PfoToPfoMap &childToParentPfoMap,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &childPfo, bool &isTrueLink, 
        bool &trueChildOrientation);

    void FillHierarchyMap(const pandora::Algorithm *const pAlgorithm, PfoToMCParticleMap &pfoToMCParticleMap,
        PfoToPfoMap &childToParentPfoMap) const;

private:

    typedef std::map<const pandora::MCParticle*, const pandora::MCParticle*> MCToMCMap;

    bool IsUpstreamTrueVertex(const PfoToMCParticleMap &pfoToMCParticleMap, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CartesianVector &upstreamVertex, const pandora::CartesianVector &downstreamVertex);

    bool GetMCParticleList(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticleList *& pMCParticleList) const;

    bool GetNeutrinoPfo(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *&pNeutrinoPfo) const;
    
    void MatchPFParticles(const pandora::Algorithm *const pAlgorithm, PfoToMCParticleMap &pfoToMCParticleMap) const;

    float GetNSpacepoints(const pandora::ParticleFlowObject *const pPfo) const;

    bool IsEMParticle(const pandora::MCParticle *const pMCParticle) const;

    const pandora::MCParticle* GetLeadEMParticle(const pandora::MCParticle *const pMCParticle) const;

    float SumEnergy(const pandora::CaloHitList &caloHitList) const;

    void GetVisibleMCHierarchy(const PfoToMCParticleMap &pfoToMCParticleMap, MCToMCMap &childToParentMCMap) const;

    void GetVisiblePfoHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const PfoToMCParticleMap &pfoToMCParticleMap,
        const MCToMCMap &childToParentMCMap, PfoToPfoMap &childToParentPfoMap) const;

    void BuildSplitHierarchy(const std::pair<const pandora::MCParticle*, pandora::PfoList> &splitPfo, PfoToPfoMap &childToParentPfoMap) const;

    float GetClosestDistance(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2) const;

    const pandora::ParticleFlowObject* BestParentInSplitHierarchy(const pandora::ParticleFlowObject *const pChildPfo, 
        const pandora::PfoList &splitParticle) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_mcParticleListName;
    std::string m_neutrinoPfoListName;
    pandora::StringVector m_pfoListNames;
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_CHEAT_HIERARCHY_TOOL_H
