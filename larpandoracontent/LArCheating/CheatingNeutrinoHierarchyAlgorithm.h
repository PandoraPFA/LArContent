/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoHierarchyAlgorithm.h
 *
 *  @brief  Header file for the cheating neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_HIERARCHY_ALGORITHM_H
#define LAR_CHEATING_NEUTRINO_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoHierarchyAlgorithm::Algorithm class
 */
class CheatingNeutrinoHierarchyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingNeutrinoHierarchyAlgorithm();

private:
    pandora::StatusCode Run();

    void GetNeutrinoPfo(const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Get the mc neutrino vector
     *
     *  @param  mcNeutrinoVector to receive the mc neutrino vector
     */
    void GetMCNeutrinoVector(pandora::MCParticleVector &mcNeutrinoVector) const;

    typedef std::unordered_map<const pandora::MCParticle *, const pandora::ParticleFlowObject *> MCParticleToPfoMap;

    /**
     *  @brief  Extract candidate daughter pfos from external lists and populate a map from main mc particle (or primary) to pfo
     *
     *  @param  mcParticleToPfoMap to receive the mc particle to pfo map
     */
    void GetMCParticleToDaughterPfoMap(MCParticleToPfoMap &mcParticleToPfoMap) const;

    /**
     *  @brief  Use information from mc particles and the mc particle to pfo map to fully-reconstruct the daughter pfo hierarchy
     *
     *  @param  pParentMCParticle the address of the (current) parent mc particle
     *  @param  pParentPfo the address of the (current) parent pfo
     *  @param  mcParticleToPfoMap the mc particle to pfo map
     */
    void CreatePfoHierarchy(const pandora::MCParticle *const pParentMCParticle, const pandora::ParticleFlowObject *const pParentPfo,
        const MCParticleToPfoMap &mcParticleToPfoMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_collapseToPrimaryMCParticles; ///< Whether to collapse mc particle hierarchies to primary particles

    std::string m_mcParticleListName;  ///< The name of the three d mc particle list name
    std::string m_neutrinoPfoListName; ///< The name of the neutrino pfo list
    pandora::StringVector m_daughterPfoListNames; ///< The list of daughter pfo list names
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_HIERARCHY_ALGORITHM_H
