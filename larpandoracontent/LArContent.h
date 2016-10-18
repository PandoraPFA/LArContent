/**
 *  @file   larpandoracontent/LArContent.h
 * 
 *  @brief  Header file detailing content for use with particle flow reconstruction at liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "larpandoracontent/LArCheating/CheatingClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArCheating/CheatingClusterCreationAlgorithm.h"
#include "larpandoracontent/LArCheating/CheatingCosmicRayIdentificationAlg.h"
#include "larpandoracontent/LArCheating/CheatingCosmicRayShowerMatchingAlg.h"
#include "larpandoracontent/LArCheating/CheatingEventSlicingTool.h"
#include "larpandoracontent/LArCheating/CheatingNeutrinoCreationAlgorithm.h"
#include "larpandoracontent/LArCheating/CheatingNeutrinoDaughterVerticesAlgorithm.h"
#include "larpandoracontent/LArCheating/CheatingPfoCharacterisationAlgorithm.h"
#include "larpandoracontent/LArCheating/CheatingPfoCreationAlgorithm.h"
#include "larpandoracontent/LArCheating/CheatingVertexCreationAlgorithm.h"

#include "larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h"
#include "larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArMonitoring/EventDisplayAlgorithm.h"
#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"
#include "larpandoracontent/LArMonitoring/ParticleAnalysisAlgorithm.h"
#include "larpandoracontent/LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.h"

#include "larpandoracontent/LArPersistency/EventReadingAlgorithm.h"
#include "larpandoracontent/LArPersistency/EventWritingAlgorithm.h"

#include "larpandoracontent/LArPlugins/LArParticleIdPlugins.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"
#include "larpandoracontent/LArStitching/StitchingObjectCreationTool.h"
#include "larpandoracontent/LArStitching/StitchingPfoMergingTool.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayTrackRecoveryAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/EventSlicingTool.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoCreationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoDaughterVerticesAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoPropertiesAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h"
#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMopUpAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerPfoMopUpAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/VertexBasedPfoMopUpAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/ParticleRecoveryAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/VertexBasedPfoRecoveryAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"
#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h"
#include "larpandoracontent/LArThreeDReco/LArShowerFragments/MopUpRemnantsTool.h"
#include "larpandoracontent/LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArShowerMatching/ClearShowersTool.h"
#include "larpandoracontent/LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.h"
#include "larpandoracontent/LArThreeDReco/LArShowerMatching/SimpleShowersTool.h"
#include "larpandoracontent/LArThreeDReco/LArShowerMatching/SplitShowersTool.h"
#include "larpandoracontent/LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.h"
//#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentGapsTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TransverseTensorVisualizationTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/UndershootTracksTool.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/SimpleClusterGrowingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/BoundedClusterMopUpAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/ConeClusterMopUpAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/IsolatedClusterMopUpAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/NearbyClusterMopUpAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/SlidingConeClusterMopUpAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/OvershootSplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArTwoDReco/TwoDParticleCreationAlgorithm.h"

#include "larpandoracontent/LArUtility/ListChangingAlgorithm.h"
#include "larpandoracontent/LArUtility/ListDissolutionAlgorithm.h"
#include "larpandoracontent/LArUtility/ListDeletionAlgorithm.h"
#include "larpandoracontent/LArUtility/ListMergingAlgorithm.h"
#include "larpandoracontent/LArUtility/ListMovingAlgorithm.h"
#include "larpandoracontent/LArUtility/ListPreparationAlgorithm.h"
#include "larpandoracontent/LArUtility/NeutrinoParentAlgorithm.h"

#include "larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.h"
#include "larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h"
#include "larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.h"

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
        d("LArCheatingPfoCharacterisation",         lar_content::CheatingPfoCharacterisationAlgorithm::Factory)                 \
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
        d("LArSlidingConePfoMopUp",                 lar_content::SlidingConePfoMopUpAlgorithm::Factory)                         \
        d("LArShowerPfoMopUp",                      lar_content::ShowerPfoMopUpAlgorithm::Factory)                              \
        d("LArVertexBasedPfoMopUp",                 lar_content::VertexBasedPfoMopUpAlgorithm::Factory)                         \
        d("LArParticleRecovery",                    lar_content::ParticleRecoveryAlgorithm::Factory)                            \
        d("LArVertexBasedPfoRecovery",              lar_content::VertexBasedPfoRecoveryAlgorithm::Factory)                      \
        d("LArThreeDRemnants",                      lar_content::ThreeDRemnantsAlgorithm::Factory)                              \
        d("LArThreeDShowers",                       lar_content::ThreeDShowersAlgorithm::Factory)                               \
        d("LArThreeDTrackFragments",                lar_content::ThreeDTrackFragmentsAlgorithm::Factory)                        \
        d("LArThreeDTransverseTracks",              lar_content::ThreeDTransverseTracksAlgorithm::Factory)                      \
        d("LArCrossGapsAssociation",                lar_content::CrossGapsAssociationAlgorithm::Factory)                        \
        d("LArCrossGapsExtension",                  lar_content::CrossGapsExtensionAlgorithm::Factory)                          \
        d("LArLongitudinalAssociation",             lar_content::LongitudinalAssociationAlgorithm::Factory)                     \
        d("LArLongitudinalExtension",               lar_content::LongitudinalExtensionAlgorithm::Factory)                       \
        d("LArSimpleClusterGrowing",                lar_content::SimpleClusterGrowingAlgorithm::Factory)                        \
        d("LArSimpleClusterMerging",                lar_content::SimpleClusterMergingAlgorithm::Factory)                        \
        d("LArTransverseAssociation",               lar_content::TransverseAssociationAlgorithm::Factory)                       \
        d("LArTransverseExtension",                 lar_content::TransverseExtensionAlgorithm::Factory)                         \
        d("LArSimpleClusterCreation",               lar_content::SimpleClusterCreationAlgorithm::Factory)                       \
        d("LArTrackClusterCreation",                lar_content::TrackClusterCreationAlgorithm::Factory)                        \
        d("LArClusteringParent",                    lar_content::ClusteringParentAlgorithm::Factory)                            \
        d("LArBoundedClusterMopUp",                 lar_content::BoundedClusterMopUpAlgorithm::Factory)                         \
        d("LArConeClusterMopUp",                    lar_content::ConeClusterMopUpAlgorithm::Factory)                            \
        d("LArIsolatedClusterMopUp",                lar_content::IsolatedClusterMopUpAlgorithm::Factory)                        \
        d("LArNearbyClusterMopUp",                  lar_content::NearbyClusterMopUpAlgorithm::Factory)                          \
        d("LArSlidingConeClusterMopUp",             lar_content::SlidingConeClusterMopUpAlgorithm::Factory)                     \
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
        d("LArTwoDParticleCreation",                lar_content::TwoDParticleCreationAlgorithm::Factory)                        \
        d("LArListChanging",                        lar_content::ListChangingAlgorithm::Factory)                                \
        d("LArListDeletion",                        lar_content::ListDeletionAlgorithm::Factory)                                \
        d("LArListDissolution",                     lar_content::ListDissolutionAlgorithm::Factory)                             \
        d("LArListMerging",                         lar_content::ListMergingAlgorithm::Factory)                                 \
        d("LArListMoving",                          lar_content::ListMovingAlgorithm::Factory)                                  \
        d("LArListPreparation",                     lar_content::ListPreparationAlgorithm::Factory)                             \
        d("LArNeutrinoParent",                      lar_content::NeutrinoParentAlgorithm::Factory)                              \
        d("LArCandidateVertexCreation",             lar_content::CandidateVertexCreationAlgorithm::Factory)                     \
        d("LArEnergyKickVertexSelection",           lar_content::EnergyKickVertexSelectionAlgorithm::Factory)                   \
        d("LArHitAngleVertexSelection",             lar_content::HitAngleVertexSelectionAlgorithm::Factory)

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
        d("LArTracksCrossingGaps",                  lar_content::TracksCrossingGapsTool::Factory)                               \
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
