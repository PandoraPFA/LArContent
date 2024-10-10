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
class VisualParticleMonitoringAlgorithm : public pandora::Algorithm
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

#ifdef MONITORING
    /**
     *  @brief  Visualize all MC particles independently.
     *
     *  This function ensures each MC particle has its own entry in the TEve tree, though common PDG colour scheme is used.
     *
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualizeIndependentMC(const LArMCParticleHelper::MCContributionMap &mcMap) const;

    /**
     *  @brief  Visualize the MC particles according to their PDG codes.
     *
     *  This function groups hits from MC particles into a single species per view.
     *
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualizeMCByPdgCode(const LArMCParticleHelper::MCContributionMap &mcMap) const;

    /**
     *  @brief  Visualize the PFO particles independently (not colour-coded by PID).
     *
     *  This function visualises each PFO.
     *
     *  @param  pfoList The list of PFOs to visualize
     **/
    void VisualizeIndependentPfo(const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Visualize the PFO particles independently (not colour-coded by PID).
     *
     *  This function visualises each PFO.
     *
     *  @param  pfoList The list of PFOs to visualize
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualizeIndependentPfo(const pandora::PfoList &pfoList, const LArMCParticleHelper::MCContributionMap &mcMap) const;

    /**
     *  @brief  Visualize the PFO particles in each slice (colour-coded by slice).
     *
     *  This function visualises each slice.
     *
     *  @param  pfoList The list of PFOs to visualize
     **/
    void VisualizeReconstructedSlice(const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Visualize the PFO particles according to their PID.
     *
     *  This function visualises each PFO according to its particle ID.
     *
     *  @param  pfoList The list of PFOs to visualize
     **/
    void VisualizePfoByParticleId(const pandora::PfoList &pfoList) const;

    /**
     *  @brief  Selects the MC particles to consider based on reconstructability criteria
     *
     *  @param  pMCList The input MC particle list
     *  @param  calotHitList The input calo hit list
     *  @param  mcMap The output map from MC particles to calo hits
     **/
    void MakeSelection(const pandora::CaloHitList *pCaloHitList,
        LArMCParticleHelper::MCContributionMap &mcMap) const;
#endif // MONITORING

    std::string m_caloHitListName;  ///< Name of input calo hit list
    std::string m_pfoListName;      ///< Name of input PFO list
    bool m_visualizeMC;             ///< Whether or not to visualize MC particles
    bool m_visualizePfo;            ///< Whether or not to visualize PFOs
    bool m_visualizeSlice;          ///< Whether or not to visualize reconstructed slices
    bool m_groupMCByPdg;            ///< Whether or not to group MC particles by particle id
    bool m_showPfoByPid;            ///< Whether or not to colour PFOs by particle id
    bool m_showPfoMatchedMC;        ///< Whether or not to display the best matched MC particle for a PFO
    bool m_isTestBeam;              ///< Whether or not this is a test beam experiment
    float m_transparencyThresholdE; ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float m_energyScaleThresholdE;  ///< Cell energy for which color is at top end of continous color palette
    float m_scalingFactor;          ///< TEve works with [cm], Pandora usually works with [mm] (but LArContent went with cm too)
};

} // namespace lar_content

#endif // LAR_VISUAL_PARTICLE_MONITORING_ALGORITHM_H
