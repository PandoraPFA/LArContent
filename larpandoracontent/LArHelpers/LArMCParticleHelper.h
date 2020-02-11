/**
 *  @file   larpandoracontent/LArHelpers/LArMCParticleHelper.h
 *
 *  @brief  Header file for the lar monte carlo particle helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_PARTICLE_HELPER_H
#define LAR_MC_PARTICLE_HELPER_H 1
#include "Pandora/PandoraInternal.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include <unordered_map>
#include <functional>

namespace lar_content
{

/**
 *  @brief  LArMCParticleHelper class
 */
class LArMCParticleHelper
{
public:
    typedef std::unordered_map<const pandora::MCParticle*, const pandora::MCParticle*> MCRelationMap;
    typedef std::unordered_map<const pandora::MCParticle*, int> MCParticleIntMap;

    typedef std::unordered_map<const pandora::MCParticle*, const pandora::ParticleFlowObject*> MCToPfoMap;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::MCParticle*> CaloHitToMCMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::ParticleFlowObject*> CaloHitToPfoMap;

    typedef std::unordered_map<const pandora::MCParticle*, pandora::CaloHitList> MCContributionMap;
    typedef std::vector<MCContributionMap> MCContributionMapVector;

    typedef std::unordered_map<const pandora::ParticleFlowObject*, pandora::CaloHitList> PfoContributionMap;
    typedef std::unordered_map<const pandora::MCParticle*, PfoContributionMap> MCToPfoMatchingMap;

    typedef std::pair<const pandora::MCParticle*, pandora::CaloHitList > MCParticleCaloHitListPair;
    typedef std::pair<const pandora::ParticleFlowObject*, pandora::CaloHitList > PfoCaloHitListPair;

    typedef std::vector<MCParticleCaloHitListPair> MCParticleToSharedHitsVector;
    typedef std::vector<PfoCaloHitListPair> PfoToSharedHitsVector;

    typedef std::map<const pandora::ParticleFlowObject*, MCParticleToSharedHitsVector> PfoToMCParticleHitSharingMap;
    typedef std::map<const pandora::MCParticle*, PfoToSharedHitsVector> MCParticleToPfoHitSharingMap;

    typedef std::pair<const pandora::ParticleFlowObject*, double> PfoCompletenessPurityPair;
    typedef std::vector<PfoCompletenessPurityPair> PfoToCompletenessPurityVector;
    typedef std::map<const pandora::MCParticle*, PfoToCompletenessPurityVector> MCParticleToPfoCompletenessPurityMap;

    /**
     *  @brief   PrimaryParameters class
     */
    class PrimaryParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        PrimaryParameters();

        unsigned int  m_minPrimaryGoodHits;       ///< the minimum number of primary good Hits
        unsigned int  m_minHitsForGoodView;       ///< the minimum number of Hits for a good view
        unsigned int  m_minPrimaryGoodViews;      ///< the minimum number of primary good views
        bool          m_selectInputHits;          ///< whether to select input hits
		bool          m_foldToPrimaries;          ///< whether to fold all hits to primary pfos and MC particles
        float         m_maxPhotonPropagation;     ///< the maximum photon propagation length
        float         m_minHitSharingFraction;    ///< the minimum Hit sharing fraction
    };

    /**
     *  @brief Fills the pfoToMCParticleCompletenessMap and pfoToMCParticlePurityMap each with vectors of pairs of MCParticles and completeness/purity
     */
    static void GetMCToPfoCompletenessPurityMaps(const MCContributionMap& mcParticleToHitsMap, const PfoContributionMap& pfoToHitsMap, const MCParticleToPfoHitSharingMap& mcParticleToPfoHitSharingMap, MCParticleToPfoCompletenessPurityMap& mcParticleToPfoCompletenessMap, MCParticleToPfoCompletenessPurityMap& mcParticleToPfoPurityMap);

    /**
     *  @brief  Returns true if passed a primary neutrino final state MCParticle
     */
    static bool IsBeamNeutrinoFinalState(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a neutrino final state MCParticle (used to prevent folding to MC primaries)
     */
    static bool IsDownstreamOfBeamNeutrino(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a primary triggered beam MCParticle
     */
    static bool IsTriggeredBeamParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a primary beam MCParticle
     */
    static bool IsBeamParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Returns true if passed a leading beam MCParticle
     */
    static bool IsLeadingBeamParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Return true if passed a primary cosmic ray MCParticle
     */
    static bool IsCosmicRay(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the nuance code of an MCParticle
     */
    static unsigned int GetNuanceCode(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is a neutrino or antineutrino
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return boolean
     */
    static bool IsNeutrino(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a provided mc particle matches the implemented definition of being primary
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return boolean
     */
    static bool IsPrimary(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a provided mc particle matches the implemented definition of being leading
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return boolean
     */
    static bool IsLeading(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Determine the position in the hierarchy for the MCParticle
     *
     *  @param  pMCParticle the address of the mc particle
     *
     *  @return integer
     */
    static int GetHierarchyTier(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Whether a mc particle is visible (i.e. long-lived charged particle)
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return boolean
     */
    static bool IsVisible(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get neutrino MC particles from an input MC particle list
     *
     *  @param  pMCParticleList the input MC particle list
     *  @param  trueNeutrinos to receive the vector of neutrino MC particles
     */
    static void GetTrueNeutrinos(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &trueNeutrinos);

    /**
     *  @brief  Get triggered test beam MC particles from an input MC particle list
     *
     *  @param  pMCParticleList the input MC particle list
     *  @param  trueTestBeamParticles to receive the vector of neutrino MC particles
     */
    static void GetTrueTestBeamParticles(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &trueTestBeamParticles);

    /**
     *  @brief  Get the primary parent mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetPrimaryMCParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get the leading particle in the hierarchy, for use at ProtoDUNE
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the primary parent mc particle
     */
    static const pandora::MCParticle *GetLeadingMCParticle(const pandora::MCParticle *const pMCParticle, const int hierarchyTierLimit = 1);

    /**
     *  @brief  Get vector of primary MC particles from an input list of MC particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcPrimaryVector the output mc particle vector
     */
    static void GetPrimaryMCParticleList(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &mcPrimaryVector);

    /**
     *  @brief  Get vector of leading MC particles from an input list of MC particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcLeadingVector the output mc particle vector
     */
    static void GetLeadingMCParticleList(const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleVector &mcLeadingVector);

    /**
     *  @brief  Get the parent mc particle
     *
     *  @param  pMCParticle the input mc particle
     *
     *  @return address of the parent mc particle
     */
    static const pandora::MCParticle *GetParentMCParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) and their primary parent mc particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcPrimaryMap the output mapping between mc particles and their parents
     */
    static void GetMCPrimaryMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) and their leading parent mc particles
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcLeadingMap the output mapping between mc particles and their leading parent
     */
    static void GetMCLeadingMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcLeadingMap);

    /**
     *  @brief  Get mapping from individual mc particles (in a provided list) to themselves
     *
     *  @param  pMCParticleList the input mc particle list
     *  @param  mcToSelfMap the output mapping between mc particles and themselves
     */
    static void GetMCToSelfMap(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcToSelfMap);

    /**
     *  @brief  Find the mc particle making the largest contribution to 2D clusters in a specified pfo
     *
     *  @param  pPfo address of the pfo to examine
     *
     *  @return address of the main mc particle
     */
    static const pandora::MCParticle *GetMainMCParticle(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Sort mc particles by their momentum
     *
     *  @param  pLhs address of first mc particle
     *  @param  pRhs address of second mc particle
     */
    static bool SortByMomentum(const pandora::MCParticle *const pLhs, const pandora::MCParticle *const pRhs);

    /**
     *  @brief  Match calo hits to their parent particles
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  mcToPrimaryMCMap input mapping between mc particles and their primaries
     *  @param  hitToPrimaryMCMap output mapping between calo hits and their main MC particle (primary MC particle if mcToPrimaryMCMap provided)
     *  @param  mcToTrueHitListMap output mapping between MC particles and their associated hits
     */
    static void GetMCParticleToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
        CaloHitToMCMap &hitToMCMap, MCContributionMap &mcToTrueHitListMap);

    /**
     *  @brief  Select primary, reconstructable mc particles that match given criteria.
     *
     *  @param  pMCParticleList the address of the list of MCParticles
     *  @param  pCaloHitList the address of the list of CaloHits
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected mcparticles to their hits
     */
    static void SelectReconstructableMCParticles(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList,
        const PrimaryParameters &parameters, std::function<bool(const pandora::MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Select leading, reconstructable mc particles in the relevant hierarchy that match given criteria.
     *
     *  @param  pMCParticleList the address of the list of MCParticles
     *  @param  pCaloHitList the address of the list of CaloHits
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected mcparticles to their hits
     */
    static void SelectReconstructableTestBeamHierarchyMCParticles(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList,
        const PrimaryParameters &parameters, std::function<bool(const pandora::MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Get unfolded mapping from Pfo to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMap the input mapping from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     */
    static void GetUnfoldedPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap, PfoContributionMap &pfoToReconstructable2DHitsMap);

    /**
     *  @brief  Get mapping from Pfo to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMap the input mapping from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     */
    static void GetPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap,
        PfoContributionMap &pfoToReconstructable2DHitsMap);

    /**
     *  @brief  Get mapping from Pfo in reconstructed test beam hierarchy to reconstructable 2D hits (=good hits belonging to a selected
     *          reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMap the input mapping from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     */
    static void GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap,
        PfoContributionMap &pfoToReconstructable2DHitsMap);

    /**
     *  @brief  Get mapping from Pfo to reconstructable 2D hits (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMaps the input vector of mappings from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     */
    static void GetPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        PfoContributionMap &pfoToReconstructable2DHitsMap);

    /**
     *  @brief  Get mapping from Pfo in reconstructed test beam hierarchy to reconstructable 2D hits (=good hits belonging to a selected
     *          reconstructable MCParticle)
     *
     *  @param  pfoList the input list of Pfos
     *  @param  selectedMCParticleToHitsMaps the input vector of mappings from selected reconstructable MCParticles to their hits
     *  @param  pfoToReconstructable2DHitsMap the output mapping from Pfos to their reconstructable 2D hits
     */
    static void GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const pandora::PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        PfoContributionMap &pfoToReconstructable2DHitsMap);

    /**
     *  @brief  Get the mappings from Pfo -> pair (reconstructable MCparticles, number of reconstructable 2D hits shared with Pfo)
     *                                reconstructable MCParticle -> pair (Pfo, number of reconstructable 2D hits shared with MCParticle)
     *
     *  @param  pfoToReconstructable2DHitsMap the input mapping from Pfos to reconstructable 2D hits
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  pfoToMCParticleHitSharingMap the output mapping from Pfos to selected reconstructable MCParticles and the number hits shared
     *  @param  mcParticleToPfoHitSharingMap the output mapping from selected reconstructable MCParticles to Pfos and the number hits shared
     */
    static void GetPfoMCParticleHitSharingMaps(const PfoContributionMap &pfoToReconstructable2DHitsMap, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        PfoToMCParticleHitSharingMap &pfoToMCParticleHitSharingMap, MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap);

    /**
     *  @brief  Select a subset of calo hits representing those that represent "reconstructable" regions of the event
     *
     *  @param  pCaloHitList the address of the input calo hit list
     *  @param  mcToPrimaryMCMap the mc particle to primary mc particle map
     *  @param  selectedCaloHitList to receive the populated selected calo hit list
     *  @param  selectInputHits whether to select input hits
     *  @param  maxPhotonPropagation the maximum photon propagation length
     */
    static void SelectCaloHits(const pandora::CaloHitList *const pCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
        pandora::CaloHitList &selectedCaloHitList, const bool selectInputHits, const float maxPhotonPropagation);

private:
    /**
     *  @brief  For a given Pfo, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pPfo the input pfo
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     */
    static void CollectReconstructable2DHits(const pandora::ParticleFlowObject *const pPfo, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        pandora::CaloHitList &reconstructableCaloHitList2D);

    /**
     *  @brief  For a given Pfo, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *          and belong in the test beam particle interaction hierarchy
     *
     *  @param  pPfo the input pfo
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     */
    static void CollectReconstructableTestBeamHierarchy2DHits(const pandora::ParticleFlowObject *const pPfo, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        pandora::CaloHitList &reconstructableCaloHitList2D);

    /**
     *  @brief  For a given Pfo list, collect the hits which are reconstructable (=good hits belonging to a selected reconstructable MCParticle)
     *
     *  @param  pfoList the input pfo list
     *  @param  selectedMCParticleToHitsMaps the input mappings from selected reconstructable MCParticles to hits
     *  @param  reconstructableCaloHitList2D the output list of reconstructable 2D calo hits in the input pfo
     */
    static void CollectReconstructable2DHits(const pandora::PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
        pandora::CaloHitList &reconstructableCaloHitList2D);

    /**
     *  @brief  Apply further selection criteria to end up with a collection of "good" calo hits that can be use to define whether
     *          a target mc particle is reconstructable.
     *
     *  @param  pSelectedCaloHitList the address of the calo hit list (typically already been through some selection procedure)
     *  @param  mcToPrimaryMCMap the mc particle to primary mc particle map
     *  @param  selectedGoodCaloHitList to receive the populated good selected calo hit list
     *  @param  selectInputHits whether to select input hits
     *  @param  minHitSharingFraction the minimum Hit sharing fraction
     */
    static void SelectGoodCaloHits(const pandora::CaloHitList *const pSelectedCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
        pandora::CaloHitList &selectedGoodCaloHitList, const bool selectInputHits, const float minHitSharingFraction);

    /**
     *  @brief  Select mc particles matching given criteria from an input list
     *
     *  @param  inputMCParticles input vector of MCParticles
     *  @param  fCriteria a function which returns a bool (= shouldSelect) for a given input MCParticle
     *  @param  selectedParticles the output vector of particles selected
     */
    static void SelectParticlesMatchingCriteria(const pandora::MCParticleVector &inputMCParticles, std::function<bool(const pandora::MCParticle *const)> fCriteria,
        pandora::MCParticleVector &selectedParticles);

    /**
     *  @brief  Filter an input vector of MCParticles to ensure they have sufficient good hits to be reconstructable
     *
     *  @param  candidateTargets candidate recontructable MCParticles
     *  @param  mcToTrueHitListMap mapping from candidates reconstructable MCParticles to their true hits
     *  @param  mcToPrimaryMCMap the mc particle to primary mc particle map
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected mcparticles to their hits
     */
    static void SelectParticlesByHitCount(const pandora::MCParticleVector &candidateTargets, const MCContributionMap &mcToTrueHitListMap,
        const MCRelationMap &mcToPrimaryMCMap, const PrimaryParameters &parameters, MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Whether it is possible to navigate from a primary mc particle to a downstream mc particle without "passing through" a neutron
     *
     *  @param  pOriginalPrimary the address of the original primary mc particle
     *  @param  pThisMCParticle the address of the current mc particle in the primary decay chain
     *  @param  pHitMCParticle the address of the mc particle associated to a calo hit
     *  @param  maxPhotonPropagation the maximum photon propagation length
     *
     *  @return boolean
     */
    static bool PassMCParticleChecks(const pandora::MCParticle *const pOriginalPrimary, const pandora::MCParticle *const pThisMCParticle,
        const pandora::MCParticle *const pHitMCParticle, const float maxPhotonPropagation);

    /**
     *  @brief  Get the hits in the intersection of two hit lists
     *
     *  @param  hitListA an input hit list
     *  @param  hitListB another input hit list
     *
     *  @return The hits that are found in both hitListA and hitListB
     */
    static pandora::CaloHitList GetSharedHits(const pandora::CaloHitList &hitListA, const pandora::CaloHitList &hitListB);
};

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_HELPER_H
