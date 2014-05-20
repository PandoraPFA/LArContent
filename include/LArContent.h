/**
 *  @file   LArContent/include/LArContent.h
 * 
 *  @brief  Header file detailing content for use with liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "LArCalculators/LArTransformationCalculator.h"
#include "LArCalculators/LArPseudoLayerCalculator.h"

#include "LArCheating/CheatingClusterCreationAlgorithm.h"
#include "LArCheating/CheatingClusterCharacterisationAlg.h"
#include "LArCheating/CheatingCosmicRayIdentificationAlg.h"
#include "LArCheating/CheatingCosmicRayShowerMatchingAlg.h"
#include "LArCheating/CheatingPfoCreationAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArMonitoring/EventDisplayAlgorithm.h"
#include "LArMonitoring/NtupleWritingAlgorithm.h"
#include "LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "LArMonitoring/VisualMonitoringAlgorithm.h"

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "LArThreeDReco/LArHitCreation/TransverseTrackHitCreationTool.h"
#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"
#include "LArThreeDReco/LArTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"
#include "LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.h"
#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"
#include "LArThreeDReco/LArTrackMatching/ClearLongitudinalTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/ClearRemnantTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/ClearTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/LongTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/MatchedEndPointsTool.h"
#include "LArThreeDReco/LArTrackMatching/MissingTrackTool.h"
#include "LArThreeDReco/LArTrackMatching/MissingTrackSegmentTool.h"
#include "LArThreeDReco/LArTrackMatching/OvershootTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/TensorVisualizationTool.h"
#include "LArThreeDReco/LArTrackMatching/TrackSplittingTool.h"
#include "LArThreeDReco/LArTrackMatching/UndershootTracksTool.h"

#include "LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/ClusterCreationAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerMergingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/TrackCleaningAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"
#include "LArTwoDReco/LArSeedFinding/SeedBranchGrowingAlgorithm.h"
#include "LArTwoDReco/LArSeedFinding/VertexSeedFindingAlgorithm.h"
#include "LArTwoDReco/TwoDParticleCreationAlgorithm.h"

#include "LArUtility/ListChangingAlgorithm.h"
#include "LArUtility/ListDissolutionAlgorithm.h"
#include "LArUtility/ListMergingAlgorithm.h"
#include "LArUtility/ListPreparationAlgorithm.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArEventDisplay",                        lar::EventDisplayAlgorithm::Factory)                                        \
        d("LArNtupleWriting",                       lar::NtupleWritingAlgorithm::Factory)                                       \
        d("LArParticleMonitoring",                  lar::ParticleMonitoringAlgorithm::Factory)                                  \
        d("LArVisualMonitoring",                    lar::VisualMonitoringAlgorithm::Factory)                                    \
        d("LArCheatingClusterCreation",             lar::CheatingClusterCreationAlgorithm::Factory)                             \
        d("LArCheatingClusterCharacterisation",     lar::CheatingClusterCharacterisationAlg::Factory)                           \
        d("LArCheatingCosmicRayIdentification",     lar::CheatingCosmicRayIdentificationAlg::Factory)                           \
        d("LArCheatingCosmicRayShowerMatching",     lar::CheatingCosmicRayShowerMatchingAlg::Factory)                           \
        d("LArCheatingPfoCreation",                 lar::CheatingPfoCreationAlgorithm::Factory)                                 \
        d("LArThreeDHitCreation",                   lar::ThreeDHitCreationAlgorithm::Factory)                                   \
        d("LArThreeDShowers",                       lar::ThreeDShowersAlgorithm::Factory)                                       \
        d("LArThreeDLongitudinalTracks",            lar::ThreeDLongitudinalTracksAlgorithm::Factory)                            \
        d("LArThreeDRemnantTracks",                 lar::ThreeDRemnantTracksAlgorithm::Factory)                                 \
        d("LArThreeDTransverseTracks",              lar::ThreeDTransverseTracksAlgorithm::Factory)                              \
        d("LArLongitudinalAssociation",             lar::LongitudinalAssociationAlgorithm::Factory)                             \
        d("LArLongitudinalExtension",               lar::LongitudinalExtensionAlgorithm::Factory)                               \
        d("LArTransverseAssociation",               lar::TransverseAssociationAlgorithm::Factory)                               \
        d("LArTransverseExtension",                 lar::TransverseExtensionAlgorithm::Factory)                                 \
        d("LArClusterCreation",                     lar::ClusterCreationAlgorithm::Factory)                                     \
        d("LArClusteringParent",                    lar::ClusteringParentAlgorithm::Factory)                                    \
        d("LArBoundedClusterMerging",               lar::BoundedClusterMergingAlgorithm::Factory)                               \
        d("LArConeBasedMerging",                    lar::ConeBasedMergingAlgorithm::Factory)                                    \
        d("LArIsolatedHitMerging",                  lar::IsolatedHitMergingAlgorithm::Factory)                                  \
        d("LArCosmicRayExtension",                  lar::CosmicRayExtensionAlgorithm::Factory)                                  \
        d("LArCosmicRayIdentification",             lar::CosmicRayIdentificationAlgorithm::Factory)                             \
        d("LArCosmicRayShowerMatching",             lar::CosmicRayShowerMatchingAlgorithm::Factory)                             \
        d("LArCosmicRayShowerMerging",              lar::CosmicRayShowerMergingAlgorithm::Factory)                              \
        d("LArCosmicRaySplitting",                  lar::CosmicRaySplittingAlgorithm::Factory)                                  \
        d("LArCosmicRayTrackMatching",              lar::CosmicRayTrackMatchingAlgorithm::Factory)                              \
        d("LArDeltaRayExtension",                   lar::DeltaRayExtensionAlgorithm::Factory)                                   \
        d("LArBranchSplitting",                     lar::BranchSplittingAlgorithm::Factory)                                     \
        d("LArCrossedTrackSplitting",               lar::CrossedTrackSplittingAlgorithm::Factory)                               \
        d("LArDeltaRaySplitting",                   lar::DeltaRaySplittingAlgorithm::Factory)                                   \
        d("LArKinkSplitting",                       lar::KinkSplittingAlgorithm::Factory)                                       \
        d("LArLayerSplitting",                      lar::LayerSplittingAlgorithm::Factory)                                      \
        d("LArTrackCleaning",                       lar::TrackCleaningAlgorithm::Factory)                                       \
        d("LArTrackConsolidation",                  lar::TrackConsolidationAlgorithm::Factory)                                  \
        d("LArVertexSplitting",                     lar::VertexSplittingAlgorithm::Factory)                                     \
        d("LArClusterCharacterisation",             lar::ClusterCharacterisationAlgorithm::Factory)                             \
        d("LArSeedBranchGrowing",                   lar::SeedBranchGrowingAlgorithm::Factory)                                   \
        d("LArVertexSeedFinding",                   lar::VertexSeedFindingAlgorithm::Factory)                                   \
        d("LArTwoDParticleCreationAlgorithm",       lar::TwoDParticleCreationAlgorithm::Factory)                                \
        d("LArListChanging",                        lar::ListChangingAlgorithm::Factory)                                        \
        d("LArListDissolution",                     lar::ListDissolutionAlgorithm::Factory)                                     \
        d("LArListMerging",                         lar::ListMergingAlgorithm::Factory)                                         \
        d("LArListPreparation",                     lar::ListPreparationAlgorithm::Factory)

    #define LAR_ALGORITHM_TOOL_LIST(d)                                                                                          \
        d("LArTransverseTrackHitCreation",          lar::TransverseTrackHitCreationTool::Factory)                               \
        d("LArClearLongitudinalTracks",             lar::ClearLongitudinalTracksTool::Factory)                                  \
        d("LArClearRemnantTracks",                  lar::ClearRemnantTracksTool::Factory)                                       \
        d("LArClearTracks",                         lar::ClearTracksTool::Factory)                                              \
        d("LArLongTracks",                          lar::LongTracksTool::Factory)                                               \
        d("LArMatchedEndPoints",                    lar::MatchedEndPointsTool::Factory)                                         \
        d("LArMissingTrack",                        lar::MissingTrackTool::Factory)                                             \
        d("LArMissingTrackSegment",                 lar::MissingTrackSegmentTool::Factory)                                      \
        d("LArOvershootTracks",                     lar::OvershootTracksTool::Factory)                                          \
        d("LArTensorVisualization",                 lar::TensorVisualizationTool::Factory)                                      \
        d("LArTrackSplitting",                      lar::TrackSplittingTool::Factory)                                           \
        d("LArUndershootTracks",                    lar::UndershootTracksTool::Factory)

    #define LAR_PARTICLE_ID_LIST(d)                                                                                             \
        d("LArEmShowerId",                          &lar::LArParticleIdHelper::LArEmShowerId)                                   \
        d("LArPhotonId",                            &lar::LArParticleIdHelper::LArPhotonId)                                     \
        d("LArElectronId",                          &lar::LArParticleIdHelper::LArElectronId)                                   \
        d("LArMuonId",                              &lar::LArParticleIdHelper::LArMuonId)

    #define LAR_SETTINGS_LIST(d)                                                                                                \
        d("LArClusterHelper",                       &lar::LArClusterHelper::ReadSettings)                                       \
        d("LArGeometryHelper",                      &lar::LArGeometryHelper::ReadSettings)                                      \
        d("LArParticleIdHelper",                    &lar::LArParticleIdHelper::ReadSettings)                                    \
        d("LArPointingClusterHelper",               &lar::LArPointingClusterHelper::ReadSettings)                               \
        d("LArThreeDHelper",                        &lar::LArThreeDHelper::ReadSettings)

    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(pandora::Pandora &pandora);

    /**
     *  @brief  Register all the lar content helper functions with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterHelperFunctions(pandora::Pandora &pandora);

    /**
     *  @brief  Register all the lar content functions with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterResetFunctions(pandora::Pandora &pandora);

    /**
     *  @brief  Register all the lar content functions with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  pLArPseudoLayerCalculator the address of the lar pseudo layer calculator
     */
    static pandora::StatusCode SetLArPseudoLayerCalculator(pandora::Pandora &pandora, lar::LArPseudoLayerCalculator *pLArPseudoLayerCalculator);

    /**
     *  @brief  Register all the lar content functions with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  pLArTransformationCalculator the address of the lar transformation calculator
     */
    static pandora::StatusCode SetLArTransformationCalculator(pandora::Pandora &pandora, lar::LArTransformationCalculator *pLArTransformationCalculator);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterAlgorithms(pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(PANDORA_REGISTER_ALGORITHM);
    LAR_ALGORITHM_TOOL_LIST(PANDORA_REGISTER_ALGORITHM_TOOL);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterHelperFunctions(pandora::Pandora &pandora)
{
    LAR_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);
    LAR_SETTINGS_LIST(PANDORA_REGISTER_SETTINGS);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterResetFunctions(pandora::Pandora &pandora)
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterResetFunction(pandora, &lar::LArThreeDHelper::Reset));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArPseudoLayerCalculator(pandora::Pandora &pandora, lar::LArPseudoLayerCalculator *pLArPseudoLayerCalculator)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerCalculator(pandora, pLArPseudoLayerCalculator));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, lar::LArGeometryHelper::SetLArPseudoLayerCalculator(pLArPseudoLayerCalculator));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArTransformationCalculator(pandora::Pandora &/*pandora*/, lar::LArTransformationCalculator *pLArTransformationCalculator)
{
    return lar::LArGeometryHelper::SetLArTransformationCalculator(pLArTransformationCalculator);
}

#endif // #ifndef LAR_CONTENT_H
