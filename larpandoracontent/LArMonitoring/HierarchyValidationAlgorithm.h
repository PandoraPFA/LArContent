/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h
 *
 *  @brief  Header file for the hierarchy validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_HIERARCHY_VALIDATION_ALGORITHM_H
#define LAR_HIERARCHY_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  HierarchyValidationAlgorithm class
 */
class HierarchyValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HierarchyValidationAlgorithm();

    virtual ~HierarchyValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Validate information at the level of MC nodes
     *
     *  @param  matchInfo The match info object to use for validation
     */
    void EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const;

    /**
     *  @brief  Validate information at the level of MC nodes
     *
     *  @param  matchInfo The match info object to use for validation
     */
    void MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const;

    /**
     *  @brief  Collates variables and fills ROOT tree for MC particles with matches
     *
     *  @param  matches The MCMatches object containing the matches
     *  @param  matchInfo The MatchInfo object with the full event context
     */
    void Fill(const LArHierarchyHelper::MCMatches &matches, const LArHierarchyHelper::MatchInfo &matchInfo) const;

    int m_event;                   ///< The current event
    std::string m_caloHitListName; ///< Name of input calo hit list
    std::string m_pfoListName;     ///< Name of input PFO list
    bool m_writeTree;              ///< Whether or not to output validation information to a ROOT file
    std::string m_filename;        ///< The name of the ROOT file to write
    std::string m_treename;        ///< The name of the ROOT tree to write
    bool m_foldToPrimaries;        ///< Whether or not to fold the hierarchy back to primary particles
    bool m_foldDynamic;            ///< Whether or not to fold the hierarchy dynamically
    bool m_foldToLeadingShowers;   ///< Whether or not to fold the hierarchy back to leading shower particles
    bool m_validateEvent;          ///< Whether to validate at the level of an event
    bool m_validateMC;             ///< Whether to validate at the level of MC nodes
};

} // namespace lar_content

#endif // LAR_HIERARCHY_VALIDATION_ALGORITHM_H
