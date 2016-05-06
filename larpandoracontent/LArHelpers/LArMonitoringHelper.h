/**
 *  @file   LArContent/LArHelpers/LArMonitoringHelper.h
 *
 *  @brief  Header file for the lar monitoring helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MONITORING_HELPER_H
#define LAR_MONITORING_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  LArMonitoringHelper class
 */
class LArMonitoringHelper
{
public:
    typedef std::unordered_map<const pandora::MCParticle*, const pandora::MCParticle*> MCRelationMap;

    typedef std::unordered_map<const pandora::MCParticle*, const pandora::ParticleFlowObject*> MCToPfoMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> PfoToMCMap;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::MCParticle*> CaloHitToMCMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::ParticleFlowObject*> CaloHitToPfoMap;

    typedef std::unordered_map<const pandora::MCParticle*, pandora::CaloHitList> MCContributionMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject*, pandora::CaloHitList> PfoContributionMap;

    typedef std::unordered_map<const pandora::MCParticle*, PfoContributionMap> MCToPfoMatchingMap;

    /**
     *  @brief  Extract a list of target pfos consisting of either i) primary/final-state pfos only, or ii) a full list of all
     *          non-neutrino pfos and their daughters
     *
     *  @param  inputPfoList the input pfo list
     *  @param  primaryPfosOnly whether to extract only primary Pfos - top-level Pfos and top-level daughters of top-level neutrinos
     *  @param  pfoList the receive the output pfo list
     */
    static void ExtractTargetPfos(const pandora::PfoList &inputPfoList, const bool primaryPfosOnly, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Create map from true neutrinos to reconstructed neutrinos
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  recoNeutrinos the input list of reconstructed neutrinos
     *  @param  hitToPrimaryMCMap input mapping between calo hits and their main primary MC particle
     *  @param  outputPrimaryMap ouput mapping between from true to reconstructed neutrinos
     */
    static void GetNeutrinoMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &recoNeutrinos,
        const CaloHitToMCMap &hitToPrimaryMCMap, MCToPfoMap &outputNeutrinoMap);

    /**
     *  @brief  Match calo hits to their parent particles
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  mcToPrimaryMCMap input mapping between mc particles and their primaries
     *  @param  hitToPrimaryMCMap output mapping between calo hits and their main primary MC particle
     *  @param  mcToTrueHitListMap output mapping between MC particles and their associated hits
     */
    static void GetMCParticleToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
        CaloHitToMCMap &hitToPrimaryMCMap, MCContributionMap &mcToTrueHitListMap);

    /**
     *  @brief  Match calo hits to their parent Pfos
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  pfoList the input list of Pfos
     *  @param  collapseToPrimaryPfos whether to collapse hits associated with daughter pfos back to the primary pfo
     *  @param  hitToPfoMap output mapping between calo hits and parent Pfo
     *  @param  pfoToHitListMap output mapping between Pfos and associated hits
     */
    static void GetPfoToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &pfoList,
        const bool collapseToPrimaryPfos, CaloHitToPfoMap &hitToPfoMap, PfoContributionMap &pfoToHitListMap);

    /**
     *  @brief  Match MC particle and Pfos
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  pfoToHitListMap the mapping between Pfos and associated hits
     *  @param  hitToPrimaryMCMap input mapping between calo hits and primary mc particles
     *  @param  mcToBestPfoMap output mapping between MC particles and best matched Pfo
     *  @param  mcToBestPfoHitsMap output mapping between MC particles and list of matched hits in best pfo
     *  @param  mcToFullPfoMatchingMap output mapping between MC particles and all matched pfos (and matched hits)
     */
    static void GetMCParticleToPfoMatches(const pandora::CaloHitList *const pCaloHitList, const PfoContributionMap &pfoToHitListMap,
        const CaloHitToMCMap &hitToPrimaryMCMap, MCToPfoMap &mcToBestPfoMap, MCContributionMap &mcToBestPfoHitsMap,
        MCToPfoMatchingMap &mcToFullPfoMatchingMap);

    /**
     *  @brief  Collect up all calo hits associated with a Pfo and its daughters
     *
     *  @param  pParentPfo the input Pfo
     *  @param  caloHitList the output calo hit list
     */
    static void CollectCaloHits(const pandora::ParticleFlowObject *const pParentPfo, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Collect up all calo hits associated with a list of Pfos and their daughters
     *
     *  @param  pfoList the input Pfo list
     *  @param  caloHitList the output calo hit list
     */
    static void CollectCaloHits(const pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Count the number of calo hits, in a provided list, of a specified type
     *
     *  @param  hitType the hit type
     *  @param  caloHitList the calo hit list
     *
     *  @return the number of calo hits of the specified type
     */
    static unsigned int CountHitsByType(const pandora::HitType hitType, const pandora::CaloHitList &caloHitList);
};

} // namespace lar_content

#endif // #ifndef LAR_MONITORING_HELPER_H
