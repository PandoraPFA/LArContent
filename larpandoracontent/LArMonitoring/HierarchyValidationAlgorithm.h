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
     *  @brief  Collates variables and fills ROOT tree for MC particles with matches
     *
     *  @param matches The MCMatches object containing the matches
     */
    void FillMatched(const LArHierarchyHelper::MCMatches &matches) const;

    /**
     *  @brief  Collates variables and fills ROOT tree for MC particles without matches
     *
     *  @param node The unmatched node
     */
    void FillUnmatchedMC(const LArHierarchyHelper::MCHierarchy::Node *pNode) const;

    /**
     *  @brief  Collates variables and fills ROOT tree for reco particles without matches
     *
     *  @param node The unmatched node
     */
    void FillUnmatchedReco(const LArHierarchyHelper::RecoHierarchy::Node *pNode) const;

    std::string m_caloHitListName; ///< Name of input calo hit list
    std::string m_pfoListName;     ///< Name of input PFO list
    bool m_writeTree;              ///< Whether or not to output validation information to a ROOT file
    std::string m_filename;        ///< The name of the ROOT file to write
    std::string m_treename;        ///< The name of the ROOT tree
    bool m_foldToPrimaries;        ///< Whether or not to fold the hierarchy back to primary particles
    bool m_foldToLeadingShowers;   ///< Whether or not to fold the hierarchy back to leading shower particles
};

} // namespace lar_content

#endif // LAR_HIERARCHY_VALIDATION_ALGORITHM_H
