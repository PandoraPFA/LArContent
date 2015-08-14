/**
 *  @file   LArContent/include/LArCheating/CheatingNeutrinoCreationAlgorithm.h
 * 
 *  @brief  Header file for the cheating neutrino creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_CREATION_ALGORITHM_H
#define LAR_CHEATING_NEUTRINO_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

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
     *  @brief  Extract candidate daughter pfos from external lists, use mc information to decided whether to add as top-level daughters of neutrino
     *          ATTN No hierarchy of the daughter pfos is currently supported (doesn't matter if collapsed to mc primaries during clustering anyway)
     * 
     *  @param  pMCNeutrino the address of the mc neutrino
     *  @param  pNeutrinoPfo the address of the neutrino pfo
     */
    void AddDaughterPfos(const pandora::MCParticle *const pMCNeutrino, const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_mcParticleListName;       ///< The name of the three d mc particle list name
    std::string             m_neutrinoPfoListName;      ///< The name of the neutrino pfo list

    std::string             m_vertexListName;           ///< The name of the neutrino vertex list
    pandora::StringVector   m_daughterPfoListNames;     ///< The list of daughter pfo list names

    float                   m_vertexTolerance;          ///< Tolerance, 3d displacement, allowed between reco neutrino vertex and mc neutrino endpoint
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingNeutrinoCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingNeutrinoCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_CREATION_ALGORITHM_H
