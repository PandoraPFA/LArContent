/**
 *  @file   larpandoracontent/LArMetrics/HierarchyValidationTool.h
 *
 *  @brief  Header file for the hierarchy validation tool class.
 *
 *  $Log: $
 */
#ifndef HIERARCHY_VALIDATION_TOOL_H
#define HIERARCHY_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{
/**
 *  @brief  HierarchyValidationTool class
 */
class HierarchyValidationTool : public BaseValidationTool
{
public:
    /**
 *  @brief  EventTreeVars struct
 */
    struct HierarchyTreeVars
    {
        int m_run;                            ///< run number
        int m_subrun;                         ///< subrun number
        int m_event;                          ///< event number
        pandora::IntVector m_trueTier;        ///< true hierarchy tier within the ‘visible’ hierarchy
        pandora::IntVector m_trueParentIndex; ///< index of the true 'visible' parent in the particle list
        pandora::IntVector m_recoTier;        ///< reco hierarchy tier within the ‘visible’ hierarchy
        pandora::IntVector m_recoParentIndex; ///< index of the reco 'visible' parent in the particle list
    };

    /**
     *  @brief  Default constructor
     */
    HierarchyValidationTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu,
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, const pandora::PfoVector &bestRecoMatch);

private:
    typedef std::map<const pandora::MCParticle *, std::pair<const pandora::MCParticle *, int>> Hierarchy;

    /**
     *  @brief  Iterative function, walks down parent-child links (skipping non-reconstructable particles)
     *          to populate the 'visible' true neutrino hierarchy
     *
     *  @param[in]   pMCParticle the current parent
     *  @param[in]   pMCParent the most recent visible ancestor
     *  @param[in]   targetMC the vector of reconstructable MCParticles
     *  @param[in]   childTier the tier of found children
     *  @param[out]  hierarchy the output hierarchy
     */
    void BuildVisibleHierarchy(const pandora::MCParticle *const pMCParticle, const pandora::MCParticle *const pMCParent,
        const pandora::MCParticleVector &targetMC, const int childTier, Hierarchy &hierarchy);

    /**
     *  @brief  Fill true hierarchy variables
     *
     *  @param[in]   pMC the MCParticle
     *  @param[in]   target the vector of reconstructable MCParticles
     *  @param[in]   hierarchy the true 'visible' hierarchy
     *  @param[out]  hierarchyTreeVars the hierarchy tree variables
     */
    void FillTrueVariables(const pandora::MCParticle *const pMC, const pandora::MCParticleVector &targetMC, const Hierarchy &hierarchy,
        HierarchyTreeVars &hierarchyTreeVars);

    /**
     *  @brief  Fill reco hierarchy variables
     *
     *  @param[in]   pBestMatch the best-matched pfo
     *  @param[in]   bestMatches the vector of best-matched pfos
     *  @param[out]  hierarchyTreeVars the hierarchy tree variables
     */
    void FillRecoVariables(const pandora::Pfo *const pBestMatch, const pandora::PfoVector &bestMatches, HierarchyTreeVars &hierarchyTreeVars);

    /**
     *  @brief  Fill the hierarchy tree
     *
     *  @param[in]  hierarchyTreeVars the hierarchy tree variables
     */
    void FillTree(HierarchyTreeVars &hierarchyTreeVars);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef HIERARCHY_VALIDATION_TOOL_H
