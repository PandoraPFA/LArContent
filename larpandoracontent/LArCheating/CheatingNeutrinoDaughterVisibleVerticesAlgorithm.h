/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoDaughterVisibleVerticesAlgorithm.h
 *
 *  @brief  Header file for the cheating neutrino daughter vertices algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_DAUGHTER_VISIBLE_VERTICES_ALGORITHM_H
#define LAR_CHEATING_NEUTRINO_DAUGHTER_VISIBLE_VERTICES_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoDaughterVisibleVerticesAlgorithm::Algorithm class
 */
class CheatingNeutrinoDaughterVisibleVerticesAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingNeutrinoDaughterVisibleVerticesAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the mapping from mc particle to primary, only required if collapsed mc particle hierarchy specified
     *
     *  @param  mcPrimaryMap to receive the mapping from mc particle to primary
     */
    void GetMCPrimaryMap(LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const;

    /**
     *  @brief  Process the list of reconstructed neutrinos
     *
     *  @param  neutrinoPfos the list of neutrino pfos
     *  @param  mcPrimaryMap the mapping from mc particle to primary, only required if collapsed mc particle hierarchy specified
     */
    void ProcessRecoNeutrinos(const pandora::PfoList &neutrinoPfos, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const;

    /**
     *  @brief  Process a daughter pfo
     *
     *  @param  pDaughterPfo the address of a daughter pfo
     *  @param  mcPrimaryMap the mapping from mc particle to primary, only required if collapsed mc particle hierarchy specified
     */
    void ProcessDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_collapseToPrimaryMCParticles; ///< Whether to collapse mc particle hierarchies to primary particles
    std::string m_mcParticleListName;    ///< The mc particle list name, required if want to collapse mc particle hierarchy

    std::string m_neutrinoListName; ///< The input list of pfo list names
    std::string m_vertexListName;   ///< The name of the output cosmic-ray vertex list
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_DAUGHTER_VERTICES_ALGORITHM_H
