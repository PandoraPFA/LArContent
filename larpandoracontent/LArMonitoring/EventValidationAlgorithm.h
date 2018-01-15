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
#include <set>
#include <vector>

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
     *  @brief  Print all the raw matching output to screen
     *
     *  @param  mcParticleToHitsMap to mc particle to hits map
     *  @param  pfoToHitsMap to pfo to hits map
     *  @param  mcParticleToPfoHitSharingMap the mc particle to pfo hit sharing map
     */
    void PrintAllOutput(const LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap, const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap,
        const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const;

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
//
//    /**
//     *  @brief  Get hits, downstream from a neutrino pfo, that are truly neutrino induced (and those that are not)
//     *
//     *  @param  pNeutrinoPfo address of the neutrino pfo
//     *  @param  neutrinoInducedHits to receive the list of neutrino-induced downstream hits
//     *  @param  otherHits to receive the list of downstream hits with non-neutrino origin
//     */
//    void GetNeutrinoHitOrigins(const pandora::Pfo *const pNeutrinoPfo, pandora::CaloHitList &neutrinoInducedHits, pandora::CaloHitList &otherHits) const;
//
//    /**
//     *  @brief  Get hits, from entire event, that are truly neutrino induced (and those that are not)
//     *
//     *  @param  neutrinoInducedHits to receive the list of neutrino-induced hits
//     *  @param  otherHits to receive the list of hits with non-neutrino origin
//     */
//    void GetEventHitOrigins(pandora::CaloHitList &neutrinoInducedHits, pandora::CaloHitList &otherHits) const;
//
//    /**
//     *  @brief  Apply a well-defined matching procedure to the comprehensive matches in the provided mc primary matching map
//     *
//     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
//     *  @param  matchingDetailsMap the matching details map, to be populated
//     */
//    void PerformMatching(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, MatchingDetailsMap &matchingDetailsMap) const;
//
//    typedef std::set<int> IntSet;
//
//    /**
//     *  @brief  Get the strongest pfo match (most matched hits) between an available mc primary and an available pfo
//     *
//     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
//     *  @param  usedMCIds the list of mc primary ids with an existing match
//     *  @param  usedPfoIds the list of pfo ids with an existing match
//     *  @param  matchingDetailsMap the matching details map, to be populated
//     */
//    bool GetStrongestPfoMatch(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, IntSet &usedMCIds, IntSet &usedPfoIds, MatchingDetailsMap &matchingDetailsMap) const;
//
//    /**
//     *  @brief  Get the best matches for any pfos left-over after the strong matching procedure
//     *
//     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
//     *  @param  usedPfoIds the list of pfo ids with an existing match
//     *  @param  matchingDetailsMap the matching details map, to be populated
//     */
//    void GetRemainingPfoMatches(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const IntSet &usedPfoIds, MatchingDetailsMap &matchingDetailsMap) const;
//
//    /**
//     *  @brief  Print the results of the matching procedure
//     *
//     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
//     *  @param  matchingDetailsMap the matching details map
//     */
//    void PrintMatchingOutput(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const;
//
//#ifdef MONITORING
//    /**
//     *  @brief  Use Pandora monitoring to visualize results of the matching procedure
//     *
//     *  @param  mcNeutrinoVector the mc neutrino vector
//     *  @param  recoNeutrinoVector the reco neutrino vector
//     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
//     *  @param  matchingDetailsMap the matching details map
//     */
//    void VisualizeMatchingOutput(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector,
//        const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const;
//
//    /**
//     *  @brief  Use Pandora monitoring to visualize vertex matches
//     *
//     *  @param  mcNeutrinoVector the mc neutrino vector
//     *  @param  recoNeutrinoVector the reco neutrino vector
//     *  @param  hitType the hitType to visualize (used for projections of 3D vertex positions)
//     */
//    void VisualizeVertexMatches(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector, const pandora::HitType hitType) const;
//
//    /**
//     *  @brief  Get name and color details for a given simple mc primary
//     *
//     *  @param  simpleMCPrimary the simple mc primary
//     *  @param  mcPrimaryMatchingMap the mc primary matching map
//     *  @param  name to receive the mc primary name
//     *  @param  colorto receive the mc primary color
//     */
//    void GetPrimaryDetails(const SimpleMCPrimary &simpleMCPrimary, const MCPrimaryMatchingMap &mcPrimaryMatchingMap, std::string &name, Color &color) const;
//
//    /**
//     *  @brief  Use Pandora monitoring to visualize particles of non-neutrino origin included as daughters of the neutrino
//     *
//     *  @param  recoNeutrinoVector the reco neutrino vector
//     *  @param  allPrimaryMatchedPfos the set of all primary matched pfos
//     *  @param  hitType the hitType to visualize
//     */
//    void VisualizeContaminants(const pandora::PfoVector &recoNeutrinoVector, const pandora::PfoSet &allPrimaryMatchedPfos, const pandora::HitType hitType) const;
//
//    /**
//     *  @brief  Use Pandora monitoring to visualize left over hits and clusters
//     *
//     *  @param  hitType the hitType to visualize, will examine all remnants of given hit type in provided hit and cluster list name(s)
//     */
//    void VisualizeRemnants(const pandora::HitType hitType) const;
//#endif
//    /**
//     *  @brief  Whether a provided mc primary passes selection, based on number of "good" hits
//     *
//     *  @param  simpleMCPrimary the simple mc primary
//     *
//     *  @return boolean
//     */
//    bool IsGoodMCPrimary(const SimpleMCPrimary &simpleMCPrimary) const;
//
//    /**
//     *  @brief  Whether a provided mc primary has a match, of any quality (use simple matched pfo list and information in matching details map)
//     *
//     *  @param  simpleMCPrimary the simple mc primary
//     *  @param  simpleMatchedPfoList the list of simple matched pfos
//     *  @param  matchingDetailsMap the matching details map
//     *
//     *  @return boolean
//     */
//    bool HasMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfoList &simpleMatchedPfoList, const MatchingDetailsMap &matchingDetailsMap) const;
//
//    /**
//     *  @brief  Whether a provided mc primary and pfo are deemed to be a good match
//     *
//     *  @param  simpleMCPrimary the simple mc primary
//     *  @param  simpleMatchedPfo the simple matched pfo
//     *
//     *  @return boolean
//     */
//    bool IsGoodMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfo &simpleMatchedPfo) const;
//
//    /**
//     *  @brief  Get a mapping from pfo to unique (on an event-by-event basis) identifier
//     *
//     *  @param  pfoList the input pfo list
//     *  @param  pfoIdMap to receive the pfo id map
//     */
//    void GetPfoIdMap(const pandora::PfoList &pfoList, PfoIdMap &pfoIdMap) const;
//
//    /**
//     *  @brief  Sort simple mc primaries by number of mc hits
//     *
//     *  @param  lhs the left-hand side
//     *  @param  rhs the right-hand side
//     *
//     *  @return boolean
//     */
//    static bool SortSimpleMCPrimaries(const SimpleMCPrimary &lhs, const SimpleMCPrimary &rhs);
//
//    /**
//     *  @brief  Sort simple matched pfos by number of matched hits
//     *
//     *  @param  lhs the left-hand side
//     *  @param  rhs the right-hand side
//     *
//     *  @return boolean
//     */
//    static bool SortSimpleMatchedPfos(const SimpleMatchedPfo &lhs, const SimpleMatchedPfo &rhs);
//
//    /**
//     *  @brief  Sort reco neutrinos by number of hits in their largest daughter particle
//     *
//     *  @param  pLhs address the left-hand side instance
//     *  @param  pRhs address the right-hand side instance
//     *
//     *  @return boolean
//     */
//    static bool SortRecoNeutrinos(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<pandora::HitType> HitTypeVector;

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list

    pandora::StringVector   m_clusterListNames;             ///< Optional list of cluster list names to examine to find left-over, remnant clusters

//    bool                    m_integrateOverRecoNeutrinos;   ///< Whether to consider particles from all reco neutrinos
//    bool                    m_useRecoNeutrinosOnly;         ///< Whether to only consider pfos that are daughters of reco neutrinos
    bool                    m_useTrueNeutrinosOnly;         ///< Whether to consider only mc particles that were neutrino induced
//    bool                    m_primaryPfosOnly;              ///< Whether to extract only primary Pfos - top-level pfos and top-level daughters of top-level neutrinos
//    bool                    m_collapseToPrimaryPfos;        ///< Whether to collapse hits associated with daughter pfos back to the primary pfo
//
//    bool                    m_selectInputHits;              ///< Whether to use only hits passing mc-based quality (is "reconstructable") checks
//    float                   m_minHitNeutrinoWeight;         ///< Minimum fraction of energy deposited by neutrino-indiced products in a single hit
//    float                   m_minHitSharingFraction;        ///< Minimum fraction of energy deposited by selected primary in a single "good" hit
//    float                   m_maxPhotonPropagation;         ///< Maximum distance travelled by photon, downstream of a track, in mc particle hierarchy

    bool                    m_printAllToScreen;             ///< Whether to print all/raw matching details to screen
//    bool                    m_printMatchingToScreen;        ///< Whether to print matching output to screen
//
//    bool                    m_visualizeMatching;            ///< Whether to use Pandora monitoring to visualize matching output
//    bool                    m_visualizeVertices;            ///< Whether to show vertices in visualization
//    bool                    m_visualizeRemnants;            ///< Whether to show remnants in visualization
//    bool                    m_visualizeGaps;                ///< Whether to show geometry (detector gaps) in visualization

    bool                    m_writeToTree;                  ///< Whether to write all/raw matching details to tree

//    int                     m_matchingMinPrimaryHits;       ///< The minimum number of good mc primary hits used in matching scheme
//    int                     m_matchingMinHitsForGoodView;   ///< The minimum number of good mc primary hits in given view to declare view to be good
//    int                     m_matchingMinPrimaryGoodViews;  ///< The minimum number of good views for a mc primary
//
//    bool                    m_useSmallPrimaries;            ///< Whether to consider matches to mc primaries with fewer than m_matchingMinPrimaryHits
//    int                     m_matchingMinSharedHits;        ///< The minimum number of shared hits used in matching scheme
//    float                   m_matchingMinCompleteness;      ///< The minimum particle completeness to declare a match
//    float                   m_matchingMinPurity;            ///< The minimum particle purity to declare a match
//
//    float                   m_vertexVisualizationDeltaR;    ///< The vertex visualization delta r value, defining good and bad vertex matches

    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file

    int                     m_fileIdentifier;               ///< The input file identifier
    int                     m_eventNumber;                  ///< The event number
};

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_ALGORITHM_H
