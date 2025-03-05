/**
 *  @file   larpandoradlcontent/LArCheating/DLCheatHierarchyTool.h
 *
 *  @brief  Header file for the cheat hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_DL_CHEAT_HIERARCHY_TOOL_H
#define LAR_DL_CHEAT_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

using namespace lar_content;

namespace lar_dl_content
{

/**
 *   @brief  DLCheatHierarchyTool to calculate variables related to the initial shower region
 */
class DLCheatHierarchyTool : public pandora::AlgorithmTool
{
public:
    typedef std::map<const pandora::ParticleFlowObject *, const pandora::MCParticle *> PfoToMCParticleMap;
    typedef std::map<const pandora::ParticleFlowObject *, std::pair<const pandora::ParticleFlowObject *, int>> ChildToParentPfoMap;

    /**
     *  @brief  Default constructor
     */
    DLCheatHierarchyTool();

    pandora::StatusCode Run(const PfoToMCParticleMap &pfoToMCParticleMap, const ChildToParentPfoMap &childToParentPfoMap,
        const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo, bool &isTrueLink, bool &trueParentOrientation, bool &trueChildOrientation);

    pandora::StatusCode Run(const PfoToMCParticleMap &pfoToMCParticleMap, const ChildToParentPfoMap &childToParentPfoMap,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &childPfo, bool &isTrueLink, bool &trueChildOrientation);

    /**
    *  @brief  Determine the true child -> parent pfo visible matches, filling the 
    *          pfo->MCParticle matching map in the process
    *
    *  @param  pAlgorithm a pointer to the pandora algorithm
    *  @param  pfoToMCParticleMap the pfo->MCParticle matching map to fill
    *  @param  childToParentPfoMap the pfo->(parentPfo, generation) map to fill
    */
    void FillHierarchyMap(
        const pandora::Algorithm *const pAlgorithm, PfoToMCParticleMap &pfoToMCParticleMap, ChildToParentPfoMap &childToParentPfoMap) const;

    /**
    *  @brief  Whether the true invisible parent of a particle is a neutron
    *
    *  @param  pPfo a pointer to the input pfo 
    *  @param  pfoToMCParticleMap the pfo->MCParticle matching map
    *  @param  isNeutronChild the boolean to fill
    * 
    *  @return a StatusCode detailing whether an answer could be found
    */
    pandora::StatusCode IsNeutronChild(
        const pandora::ParticleFlowObject *const pPfo, const PfoToMCParticleMap &pfoToMCParticleMap, bool &isNeutronChild) const;

private:
    typedef std::map<const pandora::MCParticle *, const pandora::MCParticle *> MCToMCMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
    *  @brief  Whether the true startpoint of the particle is in the upstream 
    *          position (i.e. the endpoint closest to the neutrino vertex)
    *
    *  @param  pfoToMCParticleMap the pfo->MCParticle matching map
    *  @param  pPfo a pointer to the input pfo 
    *  @param  upstreamVertex the upstream endpoint
    *  @param  downstreamVertex the downstream endpoint
    * 
    *  @return whether the true startpoint of the particle is in the upstream position
    */
    bool IsUpstreamTrueVertex(const PfoToMCParticleMap &pfoToMCParticleMap, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CartesianVector &upstreamVertex, const pandora::CartesianVector &downstreamVertex);

    /**
    *  @brief  Get the MCParticle list from Pandora
    *
    *  @param  pAlgorithm a pointer to the pandora algorithm
    *  @param  pMCParticleList a pointer to the MCParticle list
    * 
    *  @return whether the MCParticle list could be set and is filled
    */
    bool GetMCParticleList(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticleList *&pMCParticleList) const;

    /**
    *  @brief  Get the neutrino pfo
    *
    *  @param  pAlgorithm a pointer to the pandora algorithm
    *  @param  pNeutrinoPfo a pointer to the neutrino pfo
    * 
    *  @return whether the neutrino pfo could be found
    */
    bool GetNeutrinoPfo(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
    *  @brief  Perform the pfo->MCParticle matching, folding back EM shower hierarchies
    *
    *  @param  pAlgorithm a pointer to the pandora algorithm
    *  @param  pfoToMCParticleMap the pfo->MCParticle matching map
    */
    void MatchPFParticles(const pandora::Algorithm *const pAlgorithm, PfoToMCParticleMap &pfoToMCParticleMap) const;

    /**
    *  @brief  Get the number of 3D hits owned by a pfo
    *
    *  @param  pPfo a pointer to the input pfo
    */
    float GetNSpacepoints(const pandora::ParticleFlowObject *const pPfo) const;

    /**
    *  @brief  Whether a MCParticle is an electron or photon
    *
    *  @param  pMCParticle a pointer to the input MCParticle
    */
    bool IsEMParticle(const pandora::MCParticle *const pMCParticle) const;

    /**
    *  @brief  Get the leading EM Particle in an EM MCParticle hierarchy
    *
    *  @param  pMCParticle a pointer to the input EM MCParticle
    *
    *  @return the leading EM Particle in the EM MCParticle hierarchy
    */
    const pandora::MCParticle *GetLeadEMParticle(const pandora::MCParticle *const pMCParticle) const;

    /**
    *  @brief  Sum the 'energy' of the W hits in an input CaloHitList
    *
    *  @param  caloHitList the input CaloHitList
    *
    *  @return the summed energy
    */
    float SumEnergy(const pandora::CaloHitList &caloHitList) const;

    /**
    *  @brief  Determine the true child->parent MCParticle visible matches
    *
    *  @param  pfoToMCParticleMap the pfo->MCParticle matching map
    *  @param  childToParentMCMap the true child->parent MCParticle visible match map to fill
    */
    void GetVisibleMCHierarchy(const PfoToMCParticleMap &pfoToMCParticleMap, MCToMCMap &childToParentMCMap) const;

    /**
    *  @brief  Determine the true child->parent pfo visible matches
    *
    *  @param  pNeutrinoPfo a pointer to the neutrino pfo
    *  @param  pfoToMCParticleMap the pfo->MCParticle matching map
    *  @param  childToParentMCMap the true child->parent MCParticle visible match map
    *  @param  childToParentPfoMap the true child->parent pfo visible match map to fill
    */
    void GetVisiblePfoHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const PfoToMCParticleMap &pfoToMCParticleMap,
        const MCToMCMap &childToParentMCMap, ChildToParentPfoMap &childToParentPfoMap) const;

    /**
    *  @brief  Determine the 'internal' hierachy of reconstructed particles that match 
    *          to the same MCParticle i.e. split particles
    *
    *  @param  splitPfo the (MCParticle, list of matched pfos) pair
    *  @param  childToParentPfoMap the true child->parent pfo visible match map to fill
    */
    void BuildSplitHierarchy(const std::pair<const pandora::MCParticle *, pandora::PfoList> &splitPfo, ChildToParentPfoMap &childToParentPfoMap) const;

    /**
    *  @brief  Determine the closest distance between two pfos
    *
    *  @param  pPfo1 a pointer to one pfo
    *  @param  pPfo2 a pointer to the other pfo
    * 
    *  @return the closest distance between the two pfos
    */
    float GetClosestDistance(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2) const;

    /**
    *  @brief  In cases in which the true parent pfo is ambiguous (because the parent 
    *          particle has been split), find the best parent based on proximity
    *
    *  @param  pChildPfo a pointer to the child particle
    *  @param  splitParticle the list of pfos forming the split particle
    * 
    *  @return the best parent pfo
    */
    const pandora::ParticleFlowObject *BestParentInSplitHierarchy(
        const pandora::ParticleFlowObject *const pChildPfo, const pandora::PfoList &splitParticle) const;

    /**
    *  @brief  Find the children of an input parent pfo in the ChildToParentPfoMap, 
    *          and assign their generation
    *
    *  @param  pParentPfo a pointer to the parent pfo
    *  @param  generationToFind the generation to assign
    *  @param  childToParentPfoMap the true child->parent pfo visible match map to modify
    * 
    */
    void AssignGeneration(const pandora::ParticleFlowObject *const pParentPfo, const int generationToFind, ChildToParentPfoMap &childToParentPfoMap) const;

    std::string m_mcParticleListName;     ///< the MCParticle list name
    std::string m_neutrinoPfoListName;    ///< the neutrino pfo list name
    pandora::StringVector m_pfoListNames; ///< the vector of pfo list names
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_CHEAT_HIERARCHY_TOOL_H
