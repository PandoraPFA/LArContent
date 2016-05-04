/**
 *  @file   LArContent/LArContent.h
 * 
 *  @brief  Header file detailing content for use with particle flow reconstruction at liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "larpandoracontent/LArContent/LArCheating/CheatingClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingClusterCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingCosmicRayIdentificationAlg.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingCosmicRayShowerMatchingAlg.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingEventSlicingTool.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingNeutrinoCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingNeutrinoDaughterVerticesAlgorithm.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingPfoCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArCheating/CheatingVertexCreationAlgorithm.h"

#include "larpandoracontent/LArContent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h"
#include "larpandoracontent/LArContent/LArCustomParticles/TrackParticleBuildingAlgorithm.h"

#include "larpandoracontent/LArContent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArContent/LArMonitoring/EventDisplayAlgorithm.h"
#include "larpandoracontent/LArContent/LArMonitoring/EventValidationAlgorithm.h"
#include "larpandoracontent/LArContent/LArMonitoring/ParticleAnalysisAlgorithm.h"
#include "larpandoracontent/LArContent/LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "larpandoracontent/LArContent/LArMonitoring/VisualMonitoringAlgorithm.h"

#include "larpandoracontent/LArContent/LArPersistency/EventReadingAlgorithm.h"
#include "larpandoracontent/LArContent/LArPersistency/EventWritingAlgorithm.h"

#include "larpandoracontent/LArContent/LArPlugins/LArParticleIdPlugins.h"
#include "larpandoracontent/LArContent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArContent/LArPlugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArContent/LArStitching/StitchingAlgorithm.h"
#include "larpandoracontent/LArContent/LArStitching/StitchingObjectCreationTool.h"
#include "larpandoracontent/LArContent/LArStitching/StitchingPfoMergingTool.h"

#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/CosmicRayTrackRecoveryAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/EventSlicingTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/NeutrinoCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/NeutrinoDaughterVerticesAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/NeutrinoPropertiesAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArPfoMopUp/ParticleRecoveryAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArPfoMopUp/SplitShowerMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArPfoMopUp/VertexBasedPfoRecoveryAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerFragments/MopUpRemnantsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerMatching/ClearShowersTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerMatching/SimpleShowersTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArShowerMatching/SplitShowersTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/TransverseTensorVisualizationTool.h"
#include "larpandoracontent/LArContent/LArThreeDReco/LArTransverseTrackMatching/UndershootTracksTool.h"

#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/SimpleClusterGrowingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterMopUp/ProximityBasedMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/OvershootSplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArContent/LArTwoDReco/TwoDParticleCreationAlgorithm.h"

#include "larpandoracontent/LArContent/LArUtility/ListChangingAlgorithm.h"
#include "larpandoracontent/LArContent/LArUtility/ListDissolutionAlgorithm.h"
#include "larpandoracontent/LArContent/LArUtility/ListDeletionAlgorithm.h"
#include "larpandoracontent/LArContent/LArUtility/ListMergingAlgorithm.h"
#include "larpandoracontent/LArContent/LArUtility/ListMovingAlgorithm.h"
#include "larpandoracontent/LArContent/LArUtility/ListPreparationAlgorithm.h"
#include "larpandoracontent/LArContent/LArUtility/NeutrinoParentAlgorithm.h"

#include "larpandoracontent/LArContent/LArVertex/CandidateVertexCreationAlgorithm.h"
#include "larpandoracontent/LArContent/LArVertex/VertexSelectionAlgorithm.h"

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
