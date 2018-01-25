/**
 *  @file   larpandoracontent/LArControlFlow/CosmicRayTaggingMonitoringTool.h
 *
 *  @brief  Header file for the cosmic-ray tagging monitoring tool class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TAGGING_MONITORING_TOOL_H
#define LAR_COSMIC_RAY_TAGGING_MONITORING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "PandoraMonitoringApi.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CosmicRayTaggingMonitoringTool class
 */
class CosmicRayTaggingMonitoringTool : public CosmicRayTaggingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRayTaggingMonitoringTool();

    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm);

private:
    /**
     *  @brief
     */
    enum Classification
    {
        CR_MUON,
        CR_OTHER,
        TARGET,
        FRAGMENTED,
        ABSORBED,
        MIXED,
        SPARSE,
        UNCLASSIFIED
    };
    
    typedef std::map<const pandora::ParticleFlowObject*, float > PfoToFloatMap;
    typedef std::map<const pandora::ParticleFlowObject*, Classification > PfoClassificationMap;
    
    /**
     *  @brief  Calculate metrics to classify Pfos based on the target reconstructable MCParticles with which they share hits
     *
     *  @param  hitSharingMap input mapping from Pfos to MCParticles + number of shared hits pairs
     *  @param  pfoToCaloHitListMap input mapping from Pfos to their reconstructable 2D hits
     *  @param  targetsToGoodHitsMaps input mapping from target reconstructable MCParticles (those which shouldn't be tagged) to their good hits
     *  @param  pfoSignificanceMap output mapping from Pfos to their target significance
     *  @param  pfoPurityMap output mapping from Pfos to their target purities
     *  @param  pfoClassificationMap output mapping from Pfos to their classification
     */
    void CalculatePfoMetrics(const LArMCParticleHelper::PfoToMCParticleHitSharingMap &hitSharingMap, const LArMCParticleHelper::PfoContributionMap &pfoToCaloHitListMap, const LArMCParticleHelper::MCContributionMapVector &targetsToGoodHitsMaps, PfoToFloatMap &pfoSignificanceMap, PfoToFloatMap &pfoPurityMap, PfoClassificationMap &pfoClassificationMap) const;

    /**
     *  @brief  Returns true if the main MCParticle of the supplied Pfo is a muon
     */
    bool IsMainMCParticleMuon(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Classify a pfo given some metrics
     *
     *  @param  nHits the number of reconstructable hits in the Pfo
     *  @param  significance the number of target MCParticles represented by the Pfo
     *  @param  purity the fraction of reconstructable hits in the Pfo that come from a target MCParticle
     *  @param  isMuon is the main MCParticle that the Pfo represents a muon?
     */
    Classification ClassifyPfo(const unsigned int &nHits, const float &significance, const float &purity, const bool isMuon) const;

    /**
     *  @brief  Returns a unique color for each possible Pfo classification
     */
    LArFormattingHelper::Color GetClassificationColor(const Classification &classification) const;

    /**
     *  @brief  Returns a string for each classification
     */
    std::string GetClassificationName(const Classification &classification) const;

    /**
     *  @brief  Prints a table detailing all input Pfos and their classifications
     *
     *  @param  orderedPfoVector input vector of Pfos in print order
     *  @param  pfoToReconstructable2DHitsMap input mapping from Pfos to their reconstructable 2D hits
     *  @param  pfoPurityMap input mapping from Pfos to their purity
     *  @param  pfoSignificanceMap input mapping from Pfos to their significance
     *  @param  pfoClassificationMap input mapping from Pfos to their classification
     *  @param  ambiguousPfos input list of ambiguous Pfos as (not) tagged by a previous CR tagging module
     */
    void PrintPfoTable(const pandora::PfoVector &orderedPfoVector, const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, const PfoToFloatMap &pfoPurityMap, const PfoToFloatMap &pfoSignificanceMap, const PfoClassificationMap &pfoClassificationMap, const pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Read settings
     */
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    LArMCParticleHelper::ValidationParameters  m_parameters;                 ///< Parameters used to decide when an MCParticle is reconstructable
    unsigned int                               m_minHitsToConsiderTagging;   ///< The minimum number of hits to consider a Pfo for tagging
    float                                      m_minPurity;                  ///< The minimum purity to consider a Pfo as "pure"
    float                                      m_minImpurity;                ///< The minimum impurity to consider a Pfo as "impure"
    float                                      m_minSignificance;            ///< The minimum significance to consider a Pfo as "significant"
    std::string                                m_caloHitList2D;              ///< The 2D calo hit list 
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TAGGING_MONITORING_TOOL_H
