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
class MCParticleMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCParticleMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Extract details of each mc primary in a given mc contribution map
     *
     *  @param  mcContributionMap the mc contribution map
     */
    void PrintPrimaryMCParticles(const LArMCParticleHelper::MCContributionMap &mcContributionMap) const;

    /**
     *  @brief  Print information for a given mc particle to screen
     *
     *  @param  pMCParticle the address of the mc particle
     *  @param  mcToTrueHitListMap the mc to true hit list map
     *  @param  depth the depth in the mc particle decay hierarchy
     */
    void PrintMCParticle(const pandora::MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap,
        const int depth) const;

    std::string m_caloHitListName;    ///< Name of input calo hit list
    std::string m_mcParticleListName; ///< Name of input MC particle list

    bool m_useTrueNeutrinosOnly;      ///< Whether to consider only mc particles that were neutrino induced
    unsigned int m_minHitsForDisplay; ///< Min hits associated with mc particle to warrant display to terminal
    unsigned int m_minPrimaryGoodHits;         ///< Minimum hits which are deemed useful
    unsigned int m_minHitsForGoodView;         ///< Minimum hits for a view to be deemed good 
};

} // namespace lar_content

#endif // LAR_MC_PARTICLE_MONITORING_ALGORITHM_H
