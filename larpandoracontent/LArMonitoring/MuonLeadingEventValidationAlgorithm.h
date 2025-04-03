/**
 *  @file   larpandoracontent/LArMonitoring/MuonLeadingEventValidationAlgorithm.h
 *
 *  @brief  Header file for the muon leading event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MUON_LEADING_EVENT_VALIDATION_ALGORITHM_H
#define LAR_MUON_LEADING_EVENT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationBaseAlgorithm.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <map>

namespace lar_content
{
/**
 *  @brief  MuonLeadingEventValidationAlgorithm class
 */
class MuonLeadingEventValidationAlgorithm : public EventValidationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MuonLeadingEventValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~MuonLeadingEventValidationAlgorithm();

private:
    /**
     *  @brief  Fill the validation info containers
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  validationInfo to receive the validation info
     */
    void FillValidationInfo(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, ValidationInfo &validationInfo) const;

    /**
     *  @brief  Determine all reconstructable hits in cosmic ray pfos
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  recoCosmicRayHitList the output list of cosmic ray pfo reconstructable hits
     */
    void GetRecoCosmicRayHits(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, pandora::CaloHitList &recoCosmicRayHitList) const;

    /**
     *  @brief  Perform the main matching procedure
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  recoCosmicRayHitList the list of cosmic ray pfo reconstructable hits to remove from leading particle maps
     *  @param  minHitSharingFraction the minimum hit share fraction of a reconstructable hit
     *  @param  validationInfo to receive the validation info
     */
    void PerformUnfoldedMatching(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, const pandora::CaloHitList &recoCosmicRayHitList, const float minHitSharingFraction,
        ValidationInfo &validationInfo) const;

    /**
     *  @brief  Remove incorrectly reconstructed cosmic rays from main matching maps
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  validationInfo to receive the updated validation info
     */
    void RemoveIncorrectlyReconstructedCosmicRays(const pandora::MCParticleList *const pMCParticleList,
        const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList *const pPfoList, ValidationInfo &validationInfo) const;

    /**
     *  @brief  Perform the cosmic ray matching procedure and identify incorrectly reconstructed cosmic rays
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  incorrectlyReconstructedCosmicRays the output list of incorrectly reconstructed cosmic rays
     */
    void DetermineIncorrectlyReconstructedCosmicRays(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, pandora::MCParticleList &incorrectlyReconstructedCosmicRays) const;

    /**
     *  @brief  Print matching information in a provided validation info object, and write information to tree if configured to do so
     *
     *  @param  validationInfo the validation info
     *  @param  useInterpretedMatching whether to use the interpreted (rather than raw) matching information
     *  @param  printToScreen whether to print the information to screen
     *  @param  fillTree whether to write the information to tree
     */
    void ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const;

    /**
     *  @brief  Print leading MCParticle hits
     *
     *  @param  caloHitList the list of hits to print
     *  @param  isCR whether the hits belong to a MC cosmic ray or delta ray/michel electron
     */
#ifdef MONITORING
    void PrintHits(const pandora::CaloHitList &caloHitList, const bool isCR) const;
#endif

    /**
     *  @brief  Print leading pfo hits
     *
     *  @param  totalCaloHitList the list of hits to print
     *  @param  otherShowerCaloHitList the list of hits that in truth belong to a different shower
     *  @param  otherTrackCaloHitList the list of hits that in truth belong to a cosmic ray that is not the parent
     *  @param  parentTrackCaloHitList the list of hits that in truth belong to the parent cosmic ray
     *  @param  stringTag the event display marker string
     */
#ifdef MONITORING
    void PrintHits(const pandora::CaloHitList &totalCaloHitList, const pandora::CaloHitList &otherShowerCaloHitList,
        const pandora::CaloHitList &otherTrackCaloHitList, const pandora::CaloHitList &parentTrackCaloHitList, const std::string &stringTag) const;
#endif

    /**
     *  @brief  Print hits of the parent cosmic ray
     *
     *  @param  totalCaloHitList the list of hits to print
     *  @param  leadingCaloHitList the list of hits that in truth belong to the child hierarchy
     *  @param  stringTag the event display marker string
     */
    void PrintHits(const pandora::CaloHitList &totalCaloHitList, const pandora::CaloHitList &leadingCaloHitList, const std::string &stringTag) const;

    /**
     *  @brief  Fill an input contamination hit distance vector with the closest distance of each contaminant hit to the true leading particle hits
     *
     *  @param  contaminationHits the list of contaminant hits
     *  @param  leadingMCHitList the list of true MCParticles hits
     *  @param  bestMatchContaminationHitsDistance the output contaminant hit distance vector
     */
    void FillContaminationHitsDistance(const pandora::CaloHitList &contaminationHits, const pandora::CaloHitList &leadingMCHits,
        pandora::FloatVector &bestMatchContaminationHitsDistance) const;

    /**
     *  @brief  To filter out the hits of a given type from an input list
     *
     *  @param  inputList the input list of hits
     *  @param  hitType the specified TPC view
     *  @param  outputList the output list of hits of the specified list
     */
    void GetHitsOfType(const pandora::CaloHitList &inputList, const pandora::HitType hitType, pandora::CaloHitList &outputList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    LArMuonLeadingHelper::ValidationParameters m_validationParameters; ///< The definition of a reconstructable MCParticle
    bool m_removeRecoCosmicRayHits;   ///< Whether to remove the reconstructed cosmic ray hits from leading particle metrics
    bool m_deltaRayMode;              ///< Whether to run in delta ray mode
    bool m_michelMode;                ///< Whether to run in michel mode
    int m_cosmicRaysToSkip;           ///< The number of reconstructable cosmic rays to skip
    bool m_visualize;                 ///< Whether to visualize the MC and reco leading particles
    bool m_ignoreIncorrectCosmicRays; ///< Whether to remove the leading particles with incorrrectly reconstructed parents from metrics
    bool m_writeRawMatchesToTree;     ///< Whether to write all matches to output tree
    std::vector<int> m_deltaRayIDs;   ///< If filled, to contain the list leading particles to run metrics over
};

} // namespace lar_content

#endif // LAR_MUON_LEADING_EVENT_VALIDATION_ALGORITHM_H
