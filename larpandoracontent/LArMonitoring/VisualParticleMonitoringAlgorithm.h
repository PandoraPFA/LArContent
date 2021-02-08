/**
 *  @file   larpandoracontent/LArMonitoring/VisualParticleMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_VISUAL_PARTICLE_MONITORING_ALGORITHM_H
#define LAR_VISUAL_PARTICLE_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  VisualParticleMonitoringAlgorithm class
 */
class VisualParticleMonitoringAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VisualParticleMonitoringAlgorithm();

    virtual ~VisualParticleMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualize all MC particles independently.
     *
     *  This function ensures each MC particle has its own entry in the TEve tree, though common PDG colour scheme is used.
     *
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualizeIndependentMC(const LArMCParticleHelper::MCContributionMap &mcMap);

    /**
     *  @brief  Visualize the MC particles according to their PDG codes.
     *
     *  This function groups hits from MC particles into a single species per view.
     *
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualizeMCByPdgCode(const LArMCParticleHelper::MCContributionMap &mcMap);

    /**
     *  @brief  Selects the MC particles to consider based on reconstructability criteria
     *
     *  @param  pMCList The input MC particle list
     *  @param  calotHitList The input calo hit list
     *  @param  mcMap The output map from MC particles to calo hits
     **/
    void MakeSelection(const pandora::MCParticleList *pMCList, const pandora::CaloHitList *pCaloHitList,
        LArMCParticleHelper::MCContributionMap &mcMap);

    std::string     m_caloHitListName;      ///< Name of input calo hit list
    bool            m_groupByPdg;           ///< Whether or not to group MC particles by particle id
    bool            m_isTestBeam;           ///< Whether or not this is a test beam experiment
};

} // namespace lar_content

#endif // LAR_VISUAL_PARTICLE_MONITORING_ALGORITHM_H

