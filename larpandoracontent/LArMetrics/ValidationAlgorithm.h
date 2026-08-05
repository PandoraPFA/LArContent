/**
 *  @file   larpandoracontent/LArCheating/ValidationAlgorithm.h
 *
 *  @brief  Header file for the validation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VALIDATION_ALGORITHM_H
#define LAR_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  ValidationAlgorithm class
 */
class ValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ValidationAlgorithm();

    ~ValidationAlgorithm();

private:
    typedef std::vector<BaseValidationTool *> ValidationToolVector;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Whether an MCNode corresponds to a reconstructable MCParticle
     *
     *  @param[in]  pMCNode the MCNode
     *
     *  @return whether the MCParticle is deemed to be reconstructable
     */
    bool IsReconstructable(const LArHierarchyHelper::MCHierarchy::Node *pMCNode);

    std::string m_caloHitListName;               ///< The calo hit list name
    std::string m_mcParticleListName;            ///< The MCParticle list name
    pandora::StringVector m_pfoListNames;        ///< The pfo list names
    std::string m_fileName;                      ///< The name of the output validation file
    float m_minPurity;                           ///< Minimum purity to tag a node as being of good quality
    float m_minCompleteness;                     ///< Minimum completeness to tag a node as being of good quality
    unsigned int m_minRecoHits;                  ///< Minimum number of reconstructed primary good hits
    unsigned int m_minRecoHitsPerView;           ///< Minimum number of reconstructed hits for a good view
    unsigned int m_minRecoGoodViews;             ///< Minimum number of reconstructed primary good views
    bool m_removeRecoNeutrons;                   ///< Whether to remove reconstructed neutrons and their downstream particles
    bool m_selectRecoHits;                       ///< Whether to select reco hits that overlap with the MC particle hits
    ValidationToolVector m_validationToolVector; ///< The vector of validation tools ran by the validation algorithm
};

} // namespace lar_content

#endif // #ifndef LAR_VALIDATION_ALGORITHM_H
