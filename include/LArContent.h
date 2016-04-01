/**
 *  @file   LArContent/include/LArContent.h
 * 
 *  @brief  Header file detailing content for use with particle flow reconstruction at liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "LArCheating/CheatingClusterCharacterisationAlgorithm.h"
#include "LArCheating/CheatingClusterCreationAlgorithm.h"
#include "LArCheating/CheatingCosmicRayIdentificationAlg.h"
#include "LArCheating/CheatingCosmicRayShowerMatchingAlg.h"
#include "LArCheating/CheatingEventSlicingTool.h"
#include "LArCheating/CheatingNeutrinoCreationAlgorithm.h"
#include "LArCheating/CheatingNeutrinoDaughterVerticesAlgorithm.h"
#include "LArCheating/CheatingPfoCreationAlgorithm.h"
#include "LArCheating/CheatingVertexCreationAlgorithm.h"

#include "LArCustomParticles/ShowerParticleBuildingAlgorithm.h"
#include "LArCustomParticles/TrackParticleBuildingAlgorithm.h"

#include "LArHelpers/LArGeometryHelper.h"

#include "LArMonitoring/EventDisplayAlgorithm.h"
#include "LArMonitoring/EventValidationAlgorithm.h"
#include "LArMonitoring/ParticleAnalysisAlgorithm.h"
#include "LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "LArMonitoring/VisualMonitoringAlgorithm.h"

#include "LArPersistency/EventReadingAlgorithm.h"
#include "LArPersistency/EventWritingAlgorithm.h"

#include "LArPlugins/LArParticleIdPlugins.h"
#include "LArPlugins/LArPseudoLayerPlugin.h"
#include "LArPlugins/LArTransformationPlugin.h"

#include "LArStitching/StitchingAlgorithm.h"
#include "LArStitching/StitchingObjectCreationTool.h"
#include "LArStitching/StitchingPfoMergingTool.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/CosmicRayTrackRecoveryAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.h"
#include "LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"
#include "LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h"
#include "LArThreeDReco/LArEventBuilding/EventSlicingTool.h"
#include "LArThreeDReco/LArEventBuilding/NeutrinoCreationAlgorithm.h"
#include "LArThreeDReco/LArEventBuilding/NeutrinoDaughterVerticesAlgorithm.h"
#include "LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"
#include "LArThreeDReco/LArEventBuilding/NeutrinoPropertiesAlgorithm.h"
#include "LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h"
#include "LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h"
#include "LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h"
#include "LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h"
#include "LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.h"
#include "LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h"
#include "LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h"
#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h"
#include "LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"
#include "LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.h"
#include "LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h"
#include "LArThreeDReco/LArPfoMopUp/ParticleRecoveryAlgorithm.h"
#include "LArThreeDReco/LArPfoMopUp/SplitShowerMergingAlgorithm.h"
#include "LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h"
#include "LArThreeDReco/LArPfoMopUp/VertexBasedPfoRecoveryAlgorithm.h"
#include "LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"
#include "LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"
#include "LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h"
#include "LArThreeDReco/LArShowerFragments/MopUpRemnantsTool.h"
#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"
#include "LArThreeDReco/LArShowerMatching/ClearShowersTool.h"
#include "LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.h"
#include "LArThreeDReco/LArShowerMatching/SimpleShowersTool.h"
#include "LArThreeDReco/LArShowerMatching/SplitShowersTool.h"
#include "LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"
#include "LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"
#include "LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/TransverseTensorVisualizationTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/UndershootTracksTool.h"

#include "LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/SimpleClusterGrowingAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/ProximityBasedMergingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/OvershootSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"
#include "LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.h"
#include "LArTwoDReco/TwoDParticleCreationAlgorithm.h"

#include "LArUtility/ListChangingAlgorithm.h"
#include "LArUtility/ListDissolutionAlgorithm.h"
#include "LArUtility/ListDeletionAlgorithm.h"
#include "LArUtility/ListMergingAlgorithm.h"
#include "LArUtility/ListMovingAlgorithm.h"
#include "LArUtility/ListPreparationAlgorithm.h"
#include "LArUtility/NeutrinoParentAlgorithm.h"

#include "LArVertex/CandidateVertexCreationAlgorithm.h"
#include "LArVertex/VertexSelectionAlgorithm.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArEventDisplay",                        lar_content::EventDisplayAlgorithm::Factory)                                \
        d("LArEventValidation",                     lar_content::EventValidationAlgorithm::Factory)                             \
        d("LArParticleAnalysis",                    lar_content::ParticleAnalysisAlgorithm::Factory)                            \
        d("LArParticleMonitoring",                  lar_content::ParticleMonitoringAlgorithm::Factory)                          \
        d("LArVisualMonitoring",                    lar_content::VisualMonitoringAlgorithm::Factory)                            \
        d("LArEventReading",                        lar_content::EventReadingAlgorithm::Factory)                                \
        d("LArEventWriting",                        lar_content::EventWritingAlgorithm::Factory)                                \
        d("LArCheatingClusterCharacterisation",     lar_content::CheatingClusterCharacterisationAlgorithm::Factory)             \
        d("LArCheatingClusterCreation",             lar_content::CheatingClusterCreationAlgorithm::Factory)                     \
        d("LArCheatingCosmicRayIdentification",     lar_content::CheatingCosmicRayIdentificationAlg::Factory)                   \
        d("LArCheatingCosmicRayShowerMatching",     lar_content::CheatingCosmicRayShowerMatchingAlg::Factory)                   \
        d("LArCheatingNeutrinoCreation",            lar_content::CheatingNeutrinoCreationAlgorithm::Factory)                    \
        d("LArCheatingNeutrinoDaughterVertices",    lar_content::CheatingNeutrinoDaughterVerticesAlgorithm::Factory)            \
        d("LArCheatingPfoCreation",                 lar_content::CheatingPfoCreationAlgorithm::Factory)                         \
        d("LArCheatingVertexCreation",              lar_content::CheatingVertexCreationAlgorithm::Factory)                      \
        d("LArShowerParticleBuilding",              lar_content::ShowerParticleBuildingAlgorithm::Factory)                      \
        d("LArTrackParticleBuilding",               lar_content::TrackParticleBuildingAlgorithm::Factory)                       \
        d("LArStitching",                           lar_content::StitchingAlgorithm::Factory)                                   \
        d("LArCosmicRayIdentification",             lar_content::CosmicRayIdentificationAlgorithm::Factory)                     \
        d("LArCosmicRayShowerMatching",             lar_content::CosmicRayShowerMatchingAlgorithm::Factory)                     \
        d("LArCosmicRayTrackMatching",              lar_content::CosmicRayTrackMatchingAlgorithm::Factory)                      \
        d("LArCosmicRayTrackRecovery",              lar_content::CosmicRayTrackRecoveryAlgorithm::Factory)                      \
        d("LArCosmicRayVertexBuilding",             lar_content::CosmicRayVertexBuildingAlgorithm::Factory)                     \
        d("LArNeutrinoCreation",                    lar_content::NeutrinoCreationAlgorithm::Factory)                            \
        d("LArNeutrinoDaughterVertices",            lar_content::NeutrinoDaughterVerticesAlgorithm::Factory)                    \
        d("LArNeutrinoHierarchy",                   lar_content::NeutrinoHierarchyAlgorithm::Factory)                           \
        d("LArNeutrinoProperties",                  lar_content::NeutrinoPropertiesAlgorithm::Factory)                          \
        d("LArDeltaRayIdentification",              lar_content::DeltaRayIdentificationAlgorithm::Factory)                      \
        d("LArDeltaRayMatching",                    lar_content::DeltaRayMatchingAlgorithm::Factory)                            \
        d("LArThreeDHitCreation",                   lar_content::ThreeDHitCreationAlgorithm::Factory)                           \
        d("LArThreeDLongitudinalTracks",            lar_content::ThreeDLongitudinalTracksAlgorithm::Factory)                    \
        d("LArParticleRecovery",                    lar_content::ParticleRecoveryAlgorithm::Factory)                            \
        d("LArSplitShowerMerging",                  lar_content::SplitShowerMergingAlgorithm::Factory)                          \
        d("LArVertexBasedPfoMerging",               lar_content::VertexBasedPfoMergingAlgorithm::Factory)                       \
        d("LArVertexBasedPfoRecovery",              lar_content::VertexBasedPfoRecoveryAlgorithm::Factory)                      \
        d("LArThreeDRemnants",                      lar_content::ThreeDRemnantsAlgorithm::Factory)                              \
        d("LArThreeDShowers",                       lar_content::ThreeDShowersAlgorithm::Factory)                               \
        d("LArThreeDTrackFragments",                lar_content::ThreeDTrackFragmentsAlgorithm::Factory)                        \
        d("LArThreeDTransverseTracks",              lar_content::ThreeDTransverseTracksAlgorithm::Factory)                      \
        d("LArCrossGapsAssociation",                lar_content::CrossGapsAssociationAlgorithm::Factory)                        \
        d("LArLongitudinalAssociation",             lar_content::LongitudinalAssociationAlgorithm::Factory)                     \
        d("LArLongitudinalExtension",               lar_content::LongitudinalExtensionAlgorithm::Factory)                       \
        d("LArSimpleClusterGrowing",                lar_content::SimpleClusterGrowingAlgorithm::Factory)                        \
        d("LArSimpleClusterMerging",                lar_content::SimpleClusterMergingAlgorithm::Factory)                        \
        d("LArTransverseAssociation",               lar_content::TransverseAssociationAlgorithm::Factory)                       \
        d("LArTransverseExtension",                 lar_content::TransverseExtensionAlgorithm::Factory)                         \
        d("LArSimpleClusterCreation",               lar_content::SimpleClusterCreationAlgorithm::Factory)                       \
        d("LArTrackClusterCreation",                lar_content::TrackClusterCreationAlgorithm::Factory)                        \
        d("LArClusteringParent",                    lar_content::ClusteringParentAlgorithm::Factory)                            \
        d("LArBoundedClusterMerging",               lar_content::BoundedClusterMergingAlgorithm::Factory)                       \
        d("LArConeBasedMerging",                    lar_content::ConeBasedMergingAlgorithm::Factory)                            \
        d("LArIsolatedHitMerging",                  lar_content::IsolatedHitMergingAlgorithm::Factory)                          \
        d("LArProximityBasedMerging",               lar_content::ProximityBasedMergingAlgorithm::Factory)                       \
        d("LArCosmicRayExtension",                  lar_content::CosmicRayExtensionAlgorithm::Factory)                          \
        d("LArCosmicRaySplitting",                  lar_content::CosmicRaySplittingAlgorithm::Factory)                          \
        d("LArDeltaRayExtension",                   lar_content::DeltaRayExtensionAlgorithm::Factory)                           \
        d("LArDeltaRayGrowing",                     lar_content::DeltaRayGrowingAlgorithm::Factory)                             \
        d("LArBranchSplitting",                     lar_content::BranchSplittingAlgorithm::Factory)                             \
        d("LArCrossedTrackSplitting",               lar_content::CrossedTrackSplittingAlgorithm::Factory)                       \
        d("LArDeltaRaySplitting",                   lar_content::DeltaRaySplittingAlgorithm::Factory)                           \
        d("LArKinkSplitting",                       lar_content::KinkSplittingAlgorithm::Factory)                               \
        d("LArLayerSplitting",                      lar_content::LayerSplittingAlgorithm::Factory)                              \
        d("LArOvershootSplitting",                  lar_content::OvershootSplittingAlgorithm::Factory)                          \
        d("LArTrackConsolidation",                  lar_content::TrackConsolidationAlgorithm::Factory)                          \
        d("LArVertexSplitting",                     lar_content::VertexSplittingAlgorithm::Factory)                             \
        d("LArClusterCharacterisation",             lar_content::ClusterCharacterisationAlgorithm::Factory)                     \
        d("LArShowerGrowing",                       lar_content::ShowerGrowingAlgorithm::Factory)                               \
        d("LArTwoDParticleCreationAlgorithm",       lar_content::TwoDParticleCreationAlgorithm::Factory)                        \
        d("LArListChanging",                        lar_content::ListChangingAlgorithm::Factory)                                \
        d("LArListDeletion",                        lar_content::ListDeletionAlgorithm::Factory)                                \
        d("LArListDissolution",                     lar_content::ListDissolutionAlgorithm::Factory)                             \
        d("LArListMerging",                         lar_content::ListMergingAlgorithm::Factory)                                 \
        d("LArListMoving",                          lar_content::ListMovingAlgorithm::Factory)                                  \
        d("LArListPreparation",                     lar_content::ListPreparationAlgorithm::Factory)                             \
        d("LArNeutrinoParent",                      lar_content::NeutrinoParentAlgorithm::Factory)                              \
        d("LArCandidateVertexCreation",             lar_content::CandidateVertexCreationAlgorithm::Factory)                     \
        d("LArVertexSelection",                     lar_content::VertexSelectionAlgorithm::Factory)

    #define LAR_ALGORITHM_TOOL_LIST(d)                                                                                          \
        d("LArStitchingObjectCreation",             lar_content::StitchingObjectCreationTool::Factory)                          \
        d("LArStitchingPfoMerging",                 lar_content::StitchingPfoMergingTool::Factory)                              \
        d("LArCheatingEventSlicing",                lar_content::CheatingEventSlicingTool::Factory)                             \
        d("LArBranchAssociatedPfos",                lar_content::BranchAssociatedPfosTool::Factory)                             \
        d("LArEndAssociatedPfos",                   lar_content::EndAssociatedPfosTool::Factory)                                \
        d("LArEventSlicing",                        lar_content::EventSlicingTool::Factory)                                     \
        d("LArVertexAssociatedPfos",                lar_content::VertexAssociatedPfosTool::Factory)                             \
        d("LArClearShowers",                        lar_content::ClearShowersTool::Factory)                                     \
        d("LArShowerTensorVisualization",           lar_content::ShowerTensorVisualizationTool::Factory)                        \
        d("LArSimpleShowers",                       lar_content::SimpleShowersTool::Factory)                                    \
        d("LArSplitShowers",                        lar_content::SplitShowersTool::Factory)                                     \
        d("LArClearTrackFragments",                 lar_content::ClearTrackFragmentsTool::Factory)                              \
        d("LArClearLongitudinalTrackHits",          lar_content::ClearLongitudinalTrackHitsTool::Factory)                       \
        d("LArClearTransverseTrackHits",            lar_content::ClearTransverseTrackHitsTool::Factory)                         \
        d("LArDeltaRayShowerHits",                  lar_content::DeltaRayShowerHitsTool::Factory)                               \
        d("LArMultiValuedLongitudinalTrackHits",    lar_content::MultiValuedLongitudinalTrackHitsTool::Factory)                 \
        d("LArMultiValuedTransverseTrackHits",      lar_content::MultiValuedTransverseTrackHitsTool::Factory)                   \
        d("LArThreeViewShowerHits",                 lar_content::ThreeViewShowerHitsTool::Factory)                              \
        d("LArTwoViewShowerHits",                   lar_content::TwoViewShowerHitsTool::Factory)                                \
        d("LArClearLongitudinalTracks",             lar_content::ClearLongitudinalTracksTool::Factory)                          \
        d("LArMatchedEndPoints",                    lar_content::MatchedEndPointsTool::Factory)                                 \
        d("LArClearRemnants",                       lar_content::ClearRemnantsTool::Factory)                                    \
        d("LArConnectedRemnants",                   lar_content::ConnectedRemnantsTool::Factory)                                \
        d("LArMopUpRemnants",                       lar_content::MopUpRemnantsTool::Factory)                                    \
        d("LArClearTracks",                         lar_content::ClearTracksTool::Factory)                                      \
        d("LArLongTracks",                          lar_content::LongTracksTool::Factory)                                       \
        d("LArMissingTrack",                        lar_content::MissingTrackTool::Factory)                                     \
        d("LArMissingTrackSegment",                 lar_content::MissingTrackSegmentTool::Factory)                              \
        d("LArOvershootTracks",                     lar_content::OvershootTracksTool::Factory)                                  \
        d("LArTrackSplitting",                      lar_content::TrackSplittingTool::Factory)                                   \
        d("LArTransverseTensorVisualization",       lar_content::TransverseTensorVisualizationTool::Factory)                    \
        d("LArUndershootTracks",                    lar_content::UndershootTracksTool::Factory)

    #define LAR_PARTICLE_ID_LIST(d)                                                                                             \
        d("LArMuonId",                              lar_content::LArParticleIdPlugins::LArMuonId)

    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);

    /**
     *  @brief  Register the basic lar content plugins with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterBasicPlugins(const pandora::Pandora &pandora);

    /**
     *  @brief  Register pseudo layer plugin with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  pPseudoLayerPlugin the address of the pseudo layer plugin
     */
    static pandora::StatusCode SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, lar_content::LArPseudoLayerPlugin *const pLArPseudoLayerPlugin);

    /**
     *  @brief  Register lar coordinate transformation plugin with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  pLArTransformationPlugin the address of the lar transformation plugin
     */
    static pandora::StatusCode SetLArTransformationPlugin(const pandora::Pandora &pandora, lar_content::LArTransformationPlugin *const pLArTransformationPlugin);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(PANDORA_REGISTER_ALGORITHM);
    LAR_ALGORITHM_TOOL_LIST(PANDORA_REGISTER_ALGORITHM_TOOL);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterBasicPlugins(const pandora::Pandora &pandora)
{
    LAR_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, lar_content::LArPseudoLayerPlugin *const pLArPseudoLayerPlugin)
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, pLArPseudoLayerPlugin));

    return lar_content::LArGeometryHelper::SetLArPseudoLayerPlugin(pandora, pLArPseudoLayerPlugin);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArTransformationPlugin(const pandora::Pandora &pandora, lar_content::LArTransformationPlugin *const pLArTransformationPlugin)
{
    return lar_content::LArGeometryHelper::SetLArTransformationPlugin(pandora, pLArTransformationPlugin);
}

#endif // #ifndef LAR_CONTENT_H
