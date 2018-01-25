/**
 *  @file   larpandoracontent/LArMonitoring/MCParticleMonitoringAlgorithm.h
 *
 *  @brief  Header file for the mc particle monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_PARTICLE_MONITORING_ALGORITHM_H
#define LAR_MC_PARTICLE_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  MCParticleMonitoringAlgorithm class
 */
class MCParticleMonitoringAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCParticleMonitoringAlgorithm();

private:
    /**
     *  @brief SimpleMCParticle class
     */
    class SimpleMCParticle
    {
    public:
        /**
         *  @brief  Constructor
         */
        SimpleMCParticle();

        /**
         *  @brief  operator <
         * 
         *  @param  rhs object for comparison
         * 
         *  @return boolean
         */
        bool operator<(const SimpleMCParticle &rhs) const;

        int                                 m_id;                       ///< The unique identifier
        int                                 m_pdgCode;                  ///< The pdg code
        int                                 m_nMCHitsTotal;             ///< The total number of mc hits
        int                                 m_nMCHitsU;                 ///< The number of u mc hits
        int                                 m_nMCHitsV;                 ///< The number of v mc hits
        int                                 m_nMCHitsW;                 ///< The number of w mc hits
        float                               m_energy;                   ///< The energy
        pandora::CartesianVector            m_momentum;                 ///< The momentum (presumably at the vertex)
        pandora::CartesianVector            m_vertex;                   ///< The vertex
        pandora::CartesianVector            m_endpoint;                 ///< The endpoint
        const pandora::MCParticle          *m_pPandoraAddress;          ///< The address of the Pandora mc primary
    };

    typedef std::vector<SimpleMCParticle> SimpleMCParticleList;
    typedef std::unordered_map<const pandora::MCParticle*, SimpleMCParticle> SimpleMCParticleMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Extract details of each mc particle (ordered by number of true hits)
     * 
     *  @param  mcPrimaryList the mc particle list
     *  @param  mcToTrueHitListMap the mc to true hit list map
     *  @param  simpleMCParticleList to receive the populated simple mc particle list
     */
    void GetSimpleMCParticleList(const pandora::MCParticleVector &mcParticleList, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap,
        SimpleMCParticleList &simpleMCParticleList) const;

    /**
     *  @brief  Print all the mc monitoring information to screen
     * 
     *  @param  mcNeutrinoVector the mc neutrino vector
     *  @param  simpleMCParticleList the simple mc particle list
     *  @param  simpleMCParticleMap the simple mc particle map
     */
    void PrintAllOutput(const pandora::MCParticleVector &mcNeutrinoVector, const SimpleMCParticleList &simpleMCParticleList,
        const SimpleMCParticleMap &simpleMCParticleMap) const;

    /**
     *  @brief  Print information for a given mc particle to screen
     * 
     *  @param  pMCParticle the address of the mc particle
     *  @param  simpleMCParticleMap the simple mc particle map
     *  @param  depth the depth in the mc particle decay hierarchy
     */
    void PrintMCParticle(const pandora::MCParticle *const pMCParticle, const SimpleMCParticleMap &simpleMCParticleMap,
        const int depth) const;

    /**
     *  @brief  Sort simple mc particles by number of mc hits
     * 
     *  @param  lhs the left-hand side
     *  @param  rhs the right-hand side
     * 
     *  @return boolean
     */
    static bool SortSimpleMCParticles(const SimpleMCParticle &lhs, const SimpleMCParticle &rhs);

    std::string     m_caloHitListName;          ///< Name of input calo hit list
    std::string     m_mcParticleListName;       ///< Name of input MC particle list

    bool            m_useTrueNeutrinosOnly;     ///< Whether to consider only mc particles that were neutrino induced
    int             m_minHitsForDisplay;        ///< Min hits associated with mc particle to warrant display to terminal
};

} // namespace lar_content

#endif // LAR_MC_PARTICLE_MONITORING_ALGORITHM_H
