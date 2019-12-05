/**
 *  @file   larpandoracontent/LArMonitoring/ParticleEfficiencyAlgorithm.h
 *
 *  @brief Header file for the particle efficiency algorithm
 *
 * $Log: $
 */

#ifndef LAR_PARTICLE_EFFICIENCY_ALGORITHM_H
#define LAR_PARTICLE_EFFICIENCY_ALGORITHM_H

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content {

  /**
   * @brief ParticleEfficiencyAlgorithm class
   */
  class ParticleEfficiencyAlgorithm : public pandora::Algorithm {

  public:

    /**
     *  @brief Default constructor
     */
    ParticleEfficiencyAlgorithm();

    /**
     *  @brief Default destructor
     */
    ~ParticleEfficiencyAlgorithm();


    /**
     *  @brief   RecoParameters class
     */
    class RecoParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        RecoParameters();

        unsigned int  m_minPrimaryGoodHits;       ///< the minimum number of primary good Hits
        unsigned int  m_minHitsForGoodView;       ///< the minimum number of Hits for a good view
        unsigned int  m_minPrimaryGoodViews;      ///< the minimum number of primary good views
	bool          m_foldToPrimaries;          ///< whether to fold all hits to primary pfos and MC particles
        float         m_minHitSharingFraction;    ///< the minimum Hit sharing fraction
    };

  private:

    pandora::StatusCode Run();

    void FillMCToRecoHitsMap(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void SelectParticlesByHitCount(const pandora::MCParticleVector &candidateTargets, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap);

    void SelectGoodCaloHits(const pandora::CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, pandora::CaloHitList &selectedGoodCaloHitList);

    void GetMCToSelfMap(const pandora::MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToSelfMap);

    void AddMatchesEntryToTree(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void AddNoPfoEntryToTree(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_caloHitListName; // Name of input calo hit list
    std::string  m_pfoListName; // Name of input pfo list


    RecoParameters m_recoParameters;

    bool m_writeToTree;
    std::string m_treeName;
    std::string m_fileName;

    int m_eventNumber;

  };

} //namespace lar_content


#endif // LAR_PFO_VALIDATION_ALGORITHM_H
