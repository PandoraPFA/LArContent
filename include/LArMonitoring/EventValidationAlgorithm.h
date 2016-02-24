/**
 *  @file   LArContent/include/LArMonitoring/EventValidationAlgorithm.h
 *
 *  @brief  Header file for the event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_ALGORITHM_H
#define LAR_EVENT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

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

        int                     m_matchedPrimaryId;         ///< The total number of occurences
        int                     m_nMatchedHits;             ///< The number of times the primary has 0 pfo matches
        float                   m_completeness;             ///< The completeness of the match
    };

    typedef std::map<int, MatchingDetails> MatchingDetailsMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<SimpleMCPrimary, SimpleMatchedPfoList> MCPrimaryMatchingMap;
    typedef std::set<int> IntSet;

    /**
     *  @brief  Apply a well-defined matching procedure to the comprehensive matches in the provided mc primary matching map
     * 
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     */
    void PerformMatching(const MCPrimaryMatchingMap &mcPrimaryMatchingMap) const;

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

    /**
     *  @brief  Use Pandora monitoring to visualize results of the matching procedure
     * 
     *  @param  mcPrimaryMatchingMap the input/raw mc primary matching map
     *  @param  matchingDetailsMap the matching details map
     */
    void VisualizeMatchingOutput(const MCPrimaryMatchingMap &mcPrimaryMatchingMap, const MatchingDetailsMap &matchingDetailsMap) const;

    typedef std::map<const pandora::ParticleFlowObject*, int> PfoIdMap;

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

    std::string         m_caloHitListName;          ///< Name of input calo hit list
    std::string         m_mcParticleListName;       ///< Name of input MC particle list
    std::string         m_pfoListName;              ///< Name of input Pfo list

    bool                m_primaryPfosOnly;          ///< Whether to extract only primary Pfos - top-level Pfos and top-level daughters of top-level neutrinos
    bool                m_collapseToPrimaryPfos;    ///< Whether to collapse hits associated with daughter pfos back to the primary pfo

    bool                m_printAllToScreen;         ///< Whether to print all/raw matching details to screen
    bool                m_printMatchingToScreen;    ///< Whether to print matching output to screen
    bool                m_visualizeMatching;        ///< Whether to use Pandora monitoring to visualize matching output
    bool                m_writeToTree;              ///< Whether to write all/raw matching details to tree

    int                 m_matchingMinPrimaryHits;   ///< The minimum number of mc primary hits used in matching scheme
    int                 m_matchingMinSharedHits;    ///< The minimum number of shared hits used in matching scheme

    std::string         m_treeName;                 ///< Name of output tree
    std::string         m_fileName;                 ///< Name of output file

    int                 m_fileIdentifier;           ///< The input file identifier
    int                 m_eventNumber;              ///< The event number
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventValidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventValidationAlgorithm();
}

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_ALGORITHM_H
