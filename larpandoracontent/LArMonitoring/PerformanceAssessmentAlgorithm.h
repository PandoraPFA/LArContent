/**
 *  @file   larpandoracontent/LArMonitoring/PerformanceAssessmentAlgorithm.h
 *
 *  @brief Header file for the performance assessment algorithm
 *
 * $Log: $
 */

#ifndef LAR_PERFORMANCE_ASSESSMENT_ALGORITHM_H
#define LAR_PERFORMANCE_ASSESSMENT_ALGORITHM_H

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content {

  /**
   * @brief PerformanceAssessmentAlgorithm class
   */
  class PerformanceAssessmentAlgorithm : public pandora::Algorithm {

  public:
    /**
     *  @brief Default constructor
     */
    PerformanceAssessmentAlgorithm();

    ~PerformanceAssessmentAlgorithm();


  private:
    pandora::StatusCode Run();

    void FillUpstreamContainers(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleList &neutronUpstreamMCParticles, pandora::MCParticleList &photonUpstreamMCParticles, std::map<const pandora::MCParticle*, float> &photonUpstreamMCParticlesToPropagationMap);

    void FillPhotonTree(const LArMCParticleHelper::MCContributionMap &photonUpstreamMCParticleHitMap, const pandora::MCParticleList &neutronUpstreamMCParticles, const std::map<const pandora::MCParticle*, float> &photonUpstreamMCParticlesToPropagationMap);

    void FillNeutronTree(const LArMCParticleHelper::MCContributionMap &neutronUpstreamMCParticleHitMap);

    void PerformPhotonNeutronTest(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList);

    void GetMatchedMCParticlePfoCompleteness(const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const pandora::MCParticle *const pMCParticle, const pandora::ParticleFlowObject *const pPfo, double &completeness);

    void GetMatchedMCParticlePfoPurity(const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap, const pandora::MCParticle *const pMCParticle, const pandora::ParticleFlowObject *const pPfo, double &purity);

    void GetMCParticleChainVector(const pandora::MCParticle *const pMCParticle, pandora::MCParticleVector &mcParticleChainList);

    bool IsUpstream(const pandora::MCParticle *const pChainMCParticleDaughter, const pandora::MCParticle *const pEndOfChainMCParticle);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_caloHitListName; // Name of input calo hit list
    std::string  m_pfoListName; // Name of input pfo list

    bool m_writeToTree;
    bool m_printToScreen;

    std::string  m_eventTreeName;
    std::string  m_targetMCParticleTreeName;
    std::string  m_fileName;

    std::string  m_neutronFileName;
    std::string  m_photonFileName;

    int m_eventNumber;

    LArMCParticleHelper::PrimaryParameters m_parameters;

  };


} //namespace lar_content


#endif // LAR_PFO_VALIDATION_ALGORITHM_H
