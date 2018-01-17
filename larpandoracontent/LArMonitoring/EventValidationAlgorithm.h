/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationAlgorithm.h
 *
 *  @brief  Header file for the event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_ALGORITHM_H
#define LAR_EVENT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <map>

namespace lar_content
{

/**
 *  @brief  EventValidationAlgorithm class
 */
class EventValidationAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventValidationAlgorithm();

private:
    pandora::StatusCode Run();

    typedef std::unordered_map<const pandora::ParticleFlowObject*, unsigned int> PfoToIdMap;

    /**
     *  @brief  Print matching information in a provided mc particle to pfo hit sharing map
     *
     *  @param  mcParticleToHitsMap to mc particle to hits map
     *  @param  goodMCParticleToHitsMap the good mc particle to hits map
     *  @param  pfoToHitsMap to pfo to hits map
     *  @param  mcParticleToPfoHitSharingMap the mc particle to pfo hit sharing map
     *  @param  printCorrectness print number of correct elements, and corresponding fraction
     */
    void PrintOutput(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap, const LArMCParticleHelper::MCContributionMap &goodMCParticleToHitsMap,
        const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const bool printCorrectness) const;

    /**
     *  @brief  Apply an interpretative matching procedure to the comprehensive matches in the provided mc particle to pfo hit sharing map
     *
     *  @param  mcParticleToHitsMap the mc particle to hits map
     *  @param  goodMCParticleToHitsMap the good mc particle to hits map
     *  @param  pfoToHitsMap to pfo to hits map
     *  @param  mcToPfoHitSharingMap the input mc particle to pfo hit sharing map
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     */
    void InterpretMCToPfoHitSharingMap(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap, const LArMCParticleHelper::MCContributionMap &goodMCParticleToHitsMap,
        const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap,
        LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Get the strongest pfo match (most matched hits) between an available mc primary and an available pfo
     *
     *  @param  mcPrimaryVector the mc primary vector
     *  @param  mcParticleToHitsMap the mc particle to hits map
     *  @param  goodMCParticleToHitsMap the good mc particle to hits map
     *  @param  pfoToHitsMap to pfo to hits map
     *  @param  mcToPfoHitSharingMap the input mc particle to pfo hit sharing map
     *  @param  usedPfos the set of previously used pfos
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     *
     *  @return whether a strong match was identified
     */
    bool GetStrongestPfoMatch(const pandora::MCParticleVector &mcPrimaryVector, const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
        const LArMCParticleHelper::MCContributionMap &goodMCParticleToHitsMap, const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap,
        const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap, pandora::PfoSet &usedPfos,
        LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Get the best matches for any pfos left-over after the strong matching procedure
     *
     *  @param  mcPrimaryVector the mc primary vector
     *  @param  mcParticleToHitsMap the mc particle to hits map
     *  @param  goodMCParticleToHitsMap the good mc particle to hits map
     *  @param  pfoToHitsMap to pfo to hits map
     *  @param  mcToPfoHitSharingMap the input mc particle to pfo hit sharing map
     *  @param  usedPfos the set of previously used pfos
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     */
    void GetRemainingPfoMatches(const pandora::MCParticleVector &mcPrimaryVector, const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
        const LArMCParticleHelper::MCContributionMap &goodMCParticleToHitsMap, const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap,
        const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap, const pandora::PfoSet &usedPfos,
        LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Whether a provided mc primary and pfo are deemed to be a good match
     *
     *  @param  trueHits the list of true hits
     *  @param  recoHits the list of reco hits
     *  @param  sharedHits the list of shared hits
     *
     *  @return boolean
     */
    bool IsGoodMatch(const pandora::CaloHitList &trueHits, const pandora::CaloHitList &recoHits, const pandora::CaloHitList &sharedHits) const;

//    /**
//     *  @brief  Write all the raw matching output to a tree
//     *
//     *  @param  mcNeutrinoVector the mc neutrino vector
//     *  @param  recoNeutrinoVector the reco neutrino vector
//     *  @param  mcToPrimaryMCMap the mc particle to primary mc particle map
//     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
//     */
//    void WriteAllOutput(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector,
//        const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap, const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<pandora::HitType> HitTypeVector;

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list

    bool                    m_useTrueNeutrinosOnly;         ///< Whether to consider only mc particles that were neutrino induced

//    bool                    m_selectInputHits;              ///< Whether to use only hits passing mc-based quality (is "reconstructable") checks
//    float                   m_minHitNeutrinoWeight;         ///< Minimum fraction of energy deposited by neutrino-indiced products in a single hit
//    float                   m_minHitSharingFraction;        ///< Minimum fraction of energy deposited by selected primary in a single "good" hit
//    float                   m_maxPhotonPropagation;         ///< Maximum distance travelled by photon, downstream of a track, in mc particle hierarchy

    bool                    m_printAllToScreen;             ///< Whether to print all/raw matching details to screen
    bool                    m_printMatchingToScreen;        ///< Whether to print matching output to screen

    bool                    m_writeToTree;                  ///< Whether to write all/raw matching details to tree

    bool                    m_useSmallPrimaries;            ///< Whether to consider matches to mc primaries with fewer than m_matchingMinPrimaryHits
    unsigned int            m_matchingMinSharedHits;        ///< The minimum number of shared hits used in matching scheme
    float                   m_matchingMinCompleteness;      ///< The minimum particle completeness to declare a match
    float                   m_matchingMinPurity;            ///< The minimum particle purity to declare a match

    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file

    int                     m_fileIdentifier;               ///< The input file identifier
    int                     m_eventNumber;                  ///< The event number
};

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_ALGORITHM_H
