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

    void GetMatchedMCParticlePfoCompleteness(const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const pandora::MCParticle *const pMCParticle, const pandora::ParticleFlowObject *const pPfo, double &completeness);

    void GetMatchedMCParticlePfoPurity(const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap, const pandora::MCParticle *const pMCParticle, const pandora::ParticleFlowObject *const pPfo, double &purity);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_caloHitListName; // Name of input calo hit list
    std::string  m_pfoListName; // Name of input pfo list

    bool m_writeToTree;
    bool m_printToScreen;

    std::string  m_eventTreeName;
    std::string  m_targetMCParticleTreeName;
    std::string  m_fileName;

    int m_eventNumber;

    LArMCParticleHelper::PrimaryParameters m_parameters;

  };


} //namespace lar_content


#endif // LAR_PFO_VALIDATION_ALGORITHM_H
