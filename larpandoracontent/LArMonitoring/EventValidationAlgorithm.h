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

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

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
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    EventValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventValidationAlgorithm();

private:
    /**
     *  @brief SimpleMCPrimary class
     */
    class SimpleMCPrimary
    {
    public:
        /**
         *  @brief  Constructor
         */
        SimpleMCPrimary();

        /**
         *  @brief  operator <
         * 
         *  @param  rhs object for comparison
         * 
         *  @return boolean
         */
        bool operator<(const SimpleMCPrimary &rhs) const;

        int                                 m_id;                       ///< The unique identifier
        int                                 m_pdgCode;                  ///< The pdg code
        int                                 m_nMCHitsTotal;             ///< The total number of mc hits
        int                                 m_nMCHitsU;                 ///< The number of u mc hits
        int                                 m_nMCHitsV;                 ///< The number of v mc hits
        int                                 m_nMCHitsW;                 ///< The number of w mc hits
        int                                 m_nGoodMCHitsTotal;         ///< The total number of good mc hits
        int                                 m_nGoodMCHitsU;             ///< The number of good u mc hits
        int                                 m_nGoodMCHitsV;             ///< The number of good v mc hits
        int                                 m_nGoodMCHitsW;             ///< The number of good w mc hits
        float                               m_energy;                   ///< The energy
        pandora::CartesianVector            m_momentum;                 ///< The momentum (presumably at the vertex)
        pandora::CartesianVector            m_vertex;                   ///< The vertex
        pandora::CartesianVector            m_endpoint;                 ///< The endpoint
        int                                 m_nMatchedPfos;             ///< The number of matched pfos
        const pandora::MCParticle          *m_pPandoraAddress;          ///< The address of the Pandora mc primary
    };

    typedef std::vector<SimpleMCPrimary> SimpleMCPrimaryList;

    /**
     *  @brief SimpleMatchedPfo class
     */
    class SimpleMatchedPfo
    {
    public:
        /**
         *  @brief  Constructor
         */
        SimpleMatchedPfo();

        int                                 m_id;                       ///< The unique identifier
        int                                 m_parentId;                 ///< The unique identifier of the parent pfo (-1 if no parent set)
        int                                 m_pdgCode;                  ///< The pdg code
        int                                 m_nPfoHitsTotal;            ///< The total number of pfo hits
        int                                 m_nPfoHitsU;                ///< The number of u pfo hits
        int                                 m_nPfoHitsV;                ///< The number of v pfo hits
        int                                 m_nPfoHitsW;                ///< The number of w pfo hits
        int                                 m_nMatchedHitsTotal;        ///< The total number of matched hits
        int                                 m_nMatchedHitsU;            ///< The number of u matched hits
        int                                 m_nMatchedHitsV;            ///< The number of v matched hits
        int                                 m_nMatchedHitsW;            ///< The number of w matched hits
        pandora::CartesianVector            m_vertex;                   ///< The vertex (currently only filled for track pfos)
        pandora::CartesianVector            m_endpoint;                 ///< The endpoint (currently only filled for track pfos)
        pandora::CartesianVector            m_vertexDirection;          ///< The vertex direction (currently only filled for track pfos)
        pandora::CartesianVector            m_endDirection;             ///< The endpoint direction (currently only filled for track pfos)
        const pandora::ParticleFlowObject  *m_pPandoraAddress;          ///< The address of the Pandora mc primary
    };

    typedef std::vector<SimpleMatchedPfo> SimpleMatchedPfoList;

    /**
     * @brief   MatchingDetails class
     */
    class MatchingDetails
    {
    public:
        /**
         *  @brief  Default constructor
         */
        MatchingDetails();

        int                                 m_matchedPrimaryId;         ///< The total number of occurences
        int                                 m_nMatchedHits;             ///< The number of times the primary has 0 pfo matches
    };

    typedef std::map<int, MatchingDetails> MatchingDetailsMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Select a subset of calo hits representing those that represent "reconstructable" regions of the event
     * 
     *  @param  pCaloHitList the address of the input calo hit list
     *  @param  mcToPrimaryMCMap the mc particle to primary mc particle map
     *  @param  selectedCaloHitList to receive the populated selected calo hit list
     */
    void SelectCaloHits(const pandora::CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
        pandora::CaloHitList &selectedCaloHitList) const;

    /**
     *  @brief  Apply further selection criteria to end up with a collection of "good" calo hits that can be use to define whether
     *          a target mc particle is reconstructable.
     * 
     *  @param  pSelectedCaloHitList the address of the calo hit list (typically already been through some selection procedure)
     *  @param  mcToPrimaryMCMap the mc particle to primary mc particle map
     *  @param  selectedCaloHitList to receive the populated good selected calo hit list
     */
    void SelectGoodCaloHits(const pandora::CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
        pandora::CaloHitList &selectedGoodCaloHitList) const;

    /**
     *  @brief  Whether it is possible to navigate from a primary mc particle to a downstream mc particle without "passing through" a neutron
     * 
     *  @param  pOriginalPrimary the address of the original primary mc particle
     *  @param  pThisMCParticle the address of the current mc particle in the primary decay chain
     *  @param  pHitMCParticle the address of the mc particle associated to a calo hit
     * 
     *  @return boolean
     */
    bool PassMCParticleChecks(const pandora::MCParticle *const pOriginalPrimary, const pandora::MCParticle *const pThisMCParticle,
        const pandora::MCParticle *const pHitMCParticle) const;

    /**
     *  @brief  Extract details of each mc primary (ordered by number of true hits)
     * 
     *  @param  mcPrimaryList the mc primary list
     *  @param  mcToTrueHitListMap the mc to true hit list map
     *  @param  mcToGoodTrueHitListMap the mc to good true hit list map (vetoes hits with significant energy sharing)
     *  @param  mcToFullPfoMatchingMap the mc to full pfo matching map (to record number of matched pfos)
     *  @param  simpleMCPrimaryList to receive the populated simple mc primary list
     */
    void GetSimpleMCPrimaryList(const pandora::MCParticleVector &mcPrimaryList, const LArMonitoringHelper::MCContributionMap &mcToTrueHitListMap,
        const LArMonitoringHelper::MCContributionMap &mcToGoodTrueHitListMap, const LArMonitoringHelper::MCToPfoMatchingMap &mcToFullPfoMatchingMap,
        SimpleMCPrimaryList &simpleMCPrimaryList) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject*, int> PfoIdMap;
    typedef std::map<SimpleMCPrimary, SimpleMatchedPfoList> MCPrimaryMatchingMap; // SimpleMCPrimary has a defined operator<

    /**
     *  @brief  Obtain a sorted list of matched pfos for each mc primary
     * 
     *  @param  simpleMCPrimaryList the simple mc primary list
     *  @param  pfoIdMap the pfo id map
     *  @param  mcToFullPfoMatchingMap the mc to full pfo matching map
     *  @param  pfoToHitListMap the pfo to hit list map
     *  @param  mcPrimaryMatchingMap to receive the populated mc primary matching map
     */
    void GetMCPrimaryMatchingMap(const SimpleMCPrimaryList &simpleMCPrimaryList, const PfoIdMap &pfoIdMap,
        const LArMonitoringHelper::MCToPfoMatchingMap &mcToFullPfoMatchingMap, const LArMonitoringHelper::PfoContributionMap &pfoToHitListMap,
        MCPrimaryMatchingMap &mcPrimaryMatchingMap) const;

    /**
     *  @brief  Print all the raw matching output to screen
     * 
     *  @param  mcNeutrinoVector the mc neutrino vector
     *  @param  recoNeutrinoVector the reco neutrino vector
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     */
    void PrintAllOutput(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector,
        const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const;

    /**
     *  @brief  Write all the raw matching output to a tree
     * 
     *  @param  mcNeutrinoVector the mc neutrino vector
     *  @param  recoNeutrinoVector the reco neutrino vector
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     */
    void WriteAllOutput(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector,
        const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const;

    /**
     *  @brief  Apply a well-defined matching procedure to the comprehensive matches in the provided mc primary matching map
     * 
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     *  @param  matchingDetailsMap the matching details map, to be populated
     */
    void PerformMatching(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, MatchingDetailsMap &matchingDetailsMap) const;

    typedef std::set<int> IntSet;

    /**
     *  @brief  Get the strongest pfo match (most matched hits) between an available mc primary and an available pfo
     * 
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     *  @param  usedMCIds the list of mc primary ids with an existing match
     *  @param  usedPfoIds the list of pfo ids with an existing match
     *  @param  matchingDetailsMap the matching details map, to be populated
     */
    bool GetStrongestPfoMatch(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, IntSet &usedMCIds, IntSet &usedPfoIds, MatchingDetailsMap &matchingDetailsMap) const;

    /**
     *  @brief  Get the best matches for any pfos left-over after the strong matching procedure
     * 
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     *  @param  usedPfoIds the list of pfo ids with an existing match
     *  @param  matchingDetailsMap the matching details map, to be populated
     */
    void GetRemainingPfoMatches(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const IntSet &usedPfoIds, MatchingDetailsMap &matchingDetailsMap) const;

    /**
     *  @brief  Print the results of the matching procedure
     * 
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     *  @param  matchingDetailsMap the matching details map
     */
    void PrintMatchingOutput(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const;

#ifdef MONITORING
    /**
     *  @brief  Use Pandora monitoring to visualize results of the matching procedure
     * 
     *  @param  mcNeutrinoVector the mc neutrino vector
     *  @param  recoNeutrinoVector the reco neutrino vector
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     *  @param  matchingDetailsMap the matching details map
     */
    void VisualizeMatchingOutput(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector,
        const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const;

    /**
     *  @brief  Use Pandora monitoring to visualize vertex matches
     * 
     *  @param  mcNeutrinoVector the mc neutrino vector
     *  @param  recoNeutrinoVector the reco neutrino vector
     *  @param  hitType the hitType to visualize (used for projections of 3D vertex positions)
     */
    void VisualizeVertexMatches(const pandora::MCParticleVector &mcNeutrinoVector, const pandora::PfoVector &recoNeutrinoVector, const pandora::HitType hitType) const;

    /**
     *  @brief  Get name and color details for a given simple mc primary
     *
     *  @param  simpleMCPrimary the simple mc primary
     *  @param  mcPrimaryMatchingMap the mc primary matching map
     *  @param  name to receive the mc primary name
     *  @param  colorto receive the mc primary color
     */
    void GetPrimaryDetails(const SimpleMCPrimary &simpleMCPrimary, const MCPrimaryMatchingMap &mcPrimaryMatchingMap, std::string &name, Color &color) const;

    /**
     *  @brief  Use Pandora monitoring to visualize left over hits and clusters
     *
     *  @param  hitType the hitType to visualize, will examine all remnants of given hit type in provided hit and cluster list name(s)
     */
    void VisualizeRemnants(const pandora::HitType hitType) const;
#endif
    /**
     *  @brief  Whether a provided mc primary passes selection, based on number of "good" hits
     * 
     *  @param  simpleMCPrimary the simple mc primary
     * 
     *  @return boolean
     */
    bool IsGoodMCPrimary(const SimpleMCPrimary &simpleMCPrimary) const;

    /**
     *  @brief  Whether a provided mc primary has a match, of any quality (use simple matched pfo list and information in matching details map)
     * 
     *  @param  simpleMCPrimary the simple mc primary
     *  @param  simpleMatchedPfoList the list of simple matched pfos
     *  @param  matchingDetailsMap the matching details map
     * 
     *  @return boolean
     */
    bool HasMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfoList &simpleMatchedPfoList, const MatchingDetailsMap &matchingDetailsMap) const;

    /**
     *  @brief  Whether a provided mc primary and pfo are deemed to be a good match
     * 
     *  @param  simpleMCPrimary the simple mc primary
     *  @param  simpleMatchedPfo the simple matched pfo
     * 
     *  @return boolean
     */
    bool IsGoodMatch(const SimpleMCPrimary &simpleMCPrimary, const SimpleMatchedPfo &simpleMatchedPfo) const;

    /**
     *  @brief  Get a mapping from pfo to unique (on an event-by-event basis) identifier
     * 
     *  @param  pfoList the input pfo list
     *  @param  pfoIdMap to receive the pfo id map
     */
    void GetPfoIdMap(const pandora::PfoList &pfoList, PfoIdMap &pfoIdMap) const;

    /**
     *  @brief  Sort simple mc primaries by number of mc hits
     * 
     *  @param  lhs the left-hand side
     *  @param  rhs the right-hand side
     * 
     *  @return boolean
     */
    static bool SortSimpleMCPrimaries(const SimpleMCPrimary &lhs, const SimpleMCPrimary &rhs);

    /**
     *  @brief  Sort simple matched pfos by number of matched hits
     * 
     *  @param  lhs the left-hand side
     *  @param  rhs the right-hand side
     * 
     *  @return boolean
     */
    static bool SortSimpleMatchedPfos(const SimpleMatchedPfo &lhs, const SimpleMatchedPfo &rhs);

    /**
     *  @brief  Sort reco neutrinos by number of hits in their largest daughter particle
     * 
     *  @param  pLhs address the left-hand side instance
     *  @param  pRhs address the right-hand side instance
     * 
     *  @return boolean
     */
    static bool SortRecoNeutrinos(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    typedef std::vector<pandora::HitType> HitTypeVector;

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list

    pandora::StringVector   m_clusterListNames;             ///< Optional list of cluster list names to examine to find left-over, remnant clusters

    bool                    m_integrateOverSlices;          ///< Whether to consider particles from all input slices
    bool                    m_neutrinoInducedOnly;          ///< Whether to consider only mc particles that were neutrino induced
    bool                    m_primaryPfosOnly;              ///< Whether to extract only primary Pfos - top-level Pfos and top-level daughters of top-level neutrinos
    bool                    m_collapseToPrimaryPfos;        ///< Whether to collapse hits associated with daughter pfos back to the primary pfo

    float                   m_minHitSharingFraction;        ///< Minimum fraction of energy deposited by selected primary in a single "good" hit
    float                   m_maxPhotonPropagation;         ///< Maximum distance travelled by photon, downstream of a track, in mc particle hierarchy

    bool                    m_printAllToScreen;             ///< Whether to print all/raw matching details to screen
    bool                    m_printMatchingToScreen;        ///< Whether to print matching output to screen

    bool                    m_visualizeMatching;            ///< Whether to use Pandora monitoring to visualize matching output
    bool                    m_visualizeVertices;            ///< Whether to show vertices in visualization
    bool                    m_visualizeRemnants;            ///< Whether to show remnants in visualization
    bool                    m_visualizeGaps;                ///< Whether to show geometry (detector gaps) in visualization

    bool                    m_writeToTree;                  ///< Whether to write all/raw matching details to tree

    int                     m_matchingMinPrimaryHits;       ///< The minimum number of good mc primary hits used in matching scheme
    int                     m_matchingMinHitsForGoodView;   ///< The minimum number of good mc primary hits in given view to declare view to be good
    int                     m_matchingMinPrimaryGoodViews;  ///< The minimum number of good views for a mc primary

    bool                    m_useSmallPrimaries;            ///< Whether to consider matches to mc primaries with fewer than m_matchingMinPrimaryHits
    int                     m_matchingMinSharedHits;        ///< The minimum number of shared hits used in matching scheme
    float                   m_matchingMinCompleteness;      ///< The minimum particle completeness to declare a match
    float                   m_matchingMinPurity;            ///< The minimum particle purity to declare a match

    float                   m_vertexVisualizationDeltaR;    ///< The vertex visualization delta r value, defining good and bad vertex matches

    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file

    int                     m_fileIdentifier;               ///< The input file identifier
    int                     m_eventNumber;                  ///< The event number
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventValidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventValidationAlgorithm();
}

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_ALGORITHM_H
