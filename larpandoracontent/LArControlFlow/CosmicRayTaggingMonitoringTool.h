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
    LArMCParticleHelper::ValidationParameters  m_parameters;                 ///< Parameters used to decide when an MCParticle is reconstructable
    unsigned int                               m_minHitsToConsiderTagging;   ///< The minimum number of hits to consider a Pfo for tagging
    float                                      m_minPurity;                  ///< The minimum purity to consider a Pfo as "pure"
    float                                      m_minImpurity;                ///< The minimum impurity to consider a Pfo as "impure"
    float                                      m_minSignificance;            ///< The minimum significance to consider a Pfo as "significant"
    std::string                                m_caloHitList2D;              ///< The 2D calo hit list 

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
     *  @brief
     */
    void CalculatePfoMetrics(const LArMCParticleHelper::PfoToMCParticleHitSharingMap &hitSharingMap, const LArMCParticleHelper::PfoContributionMap &pfoTo2DHitsMap, const LArMCParticleHelper::MCContributionMapVector &nuMCParticlesToGoodHitsMaps, PfoToFloatMap &pfoSignificanceMap, PfoToFloatMap &pfoPurityMap, PfoClassificationMap &pfoClassificationMap) const;

    /**
     *  @brief
     */
    bool IsMainMCParticleMuon(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief
     */
    Classification ClassifyPfo(const unsigned int &nHits, const float &significance, const float &purity, const bool isMuon) const;

    /**
     *  @brief
     */
    LArFormattingHelper::Color GetClassificationColor(const Classification &classification) const;

    /**
     *  @brief
     */
    std::string GetClassificationName(const Classification &classification) const;

    /**
     *  @brief
     */
    void PrintPfoTable(const pandora::PfoVector &orderedPfoVector, const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, const PfoToFloatMap &pfoPurityMap, const PfoToFloatMap &pfoSignificanceMap, const PfoClassificationMap &pfoClassificationMap, const pandora::PfoList &ambiguousPfos) const;
    /**
     *  @brief
     */
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TAGGING_MONITORING_TOOL_H
