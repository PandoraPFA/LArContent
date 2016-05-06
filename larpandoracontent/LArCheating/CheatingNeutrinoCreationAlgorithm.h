/**
 *  @file   LArContent/LArCheating/CheatingNeutrinoCreationAlgorithm.h
 * 
 *  @brief  Header file for the cheating neutrino creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_CREATION_ALGORITHM_H
#define LAR_CHEATING_NEUTRINO_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoCreationAlgorithm::Algorithm class
 */
class CheatingNeutrinoCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    CheatingNeutrinoCreationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the address of the mc neutrino
     * 
     *  @param  pMCNeutrino to receive the address of the mc neutrino
     */
    void GetMCNeutrino(const pandora::MCParticle *&pMCNeutrino) const;

    /**
     *  @brief  Create and save a neutrino pfo with properties dictated by the mc neutrino
     * 
     *  @param  pMCNeutrino the address of the mc neutrino
     *  @param  pNeutrinoPfo to receive the address of the neutrino pfo
     */
    void CreateAndSaveNeutrinoPfo(const pandora::MCParticle *const pMCNeutrino, const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Extract reconstructed vertex from external list, check its position agrees with mc neutrino, and add to pfo
     * 
     *  @param  pMCNeutrino the address of the mc neutrino
     *  @param  pNeutrinoPfo the address of the neutrino pfo
     */
    void AddNeutrinoVertex(const pandora::MCParticle *const pMCNeutrino, const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    /**
     *  @brief  Get the mapping from mc particle to primary, only required if collapsed mc particle hierarchy specified
     * 
     *  @param  mcPrimaryMap to receive the mapping from mc particle to primary
     */
    void GetMCPrimaryMap(LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const;

    typedef std::unordered_map<const pandora::MCParticle*, const pandora::ParticleFlowObject*> MCParticleToPfoMap;

    /**
     *  @brief  Extract candidate daughter pfos from external lists and populate a map from main mc particle (or primary) to pfo
     * 
     *  @param  mcPrimaryMap the mapping from mc particle to primary, only required if collapsed mc particle hierarchy specified
     *  @param  mcParticleToPfoMap to receive the mc particle to pfo map
     */
    void GetMCParticleToDaughterPfoMap(const LArMCParticleHelper::MCRelationMap &mcPrimaryMap, MCParticleToPfoMap &mcParticleToPfoMap) const;

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

    bool                    m_collapseToPrimaryMCParticles; ///< Whether to collapse mc particle hierarchies to primary particles

    std::string             m_mcParticleListName;           ///< The name of the three d mc particle list name
    std::string             m_neutrinoPfoListName;          ///< The name of the neutrino pfo list

    std::string             m_vertexListName;               ///< The name of the neutrino vertex list
    pandora::StringVector   m_daughterPfoListNames;         ///< The list of daughter pfo list names

    float                   m_vertexTolerance;              ///< Tolerance, 3d displacement, allowed between reco neutrino vertex and mc neutrino endpoint
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingNeutrinoCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingNeutrinoCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_CREATION_ALGORITHM_H
