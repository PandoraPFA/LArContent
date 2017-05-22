/**
 *  @file   larpandoracontent/LArContent.cc
 * 
 *  @brief  Factory implementations for content intended for use with particle flow reconstruction at liquid argon time projection chambers
 * 
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"
#include "Pandora/Pandora.h"

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

#include "larpandoracontent/LArCustomParticles/PcaShowerParticleBuildingAlgorithm.h"
#include "larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h"
#include "larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"
#include "larpandoracontent/LArMonitoring/MCParticleMonitoringAlgorithm.h"
#include "larpandoracontent/LArMonitoring/VisualMonitoringAlgorithm.h"

#include "larpandoracontent/LArPersistency/EventReadingAlgorithm.h"
#include "larpandoracontent/LArPersistency/EventWritingAlgorithm.h"

#include "larpandoracontent/LArPlugins/LArParticleIdPlugins.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"
#include "larpandoracontent/LArStitching/StitchingCosmicRayMergingTool.h"
#include "larpandoracontent/LArStitching/StitchingObjectCreationTool.h"

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
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TransverseTensorVisualizationTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/UndershootTracksTool.h"

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"
#include "larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"
#include "larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.h"

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"

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
#include "larpandoracontent/LArVertex/SVMVertexSelectionAlgorithm.h"

#include "larpandoracontent/LArContent.h"

#define LAR_ALGORITHM_LIST(d)                                                                                                   \
    d("LArEventValidation",                     EventValidationAlgorithm)                                                       \
    d("LArMCParticleMonitoring",                MCParticleMonitoringAlgorithm)                                                  \
    d("LArVisualMonitoring",                    VisualMonitoringAlgorithm)                                                      \
    d("LArEventReading",                        EventReadingAlgorithm)                                                          \
    d("LArEventWriting",                        EventWritingAlgorithm)                                                          \
    d("LArCheatingClusterCharacterisation",     CheatingClusterCharacterisationAlgorithm)                                       \
    d("LArCheatingClusterCreation",             CheatingClusterCreationAlgorithm)                                               \
    d("LArCheatingCosmicRayIdentification",     CheatingCosmicRayIdentificationAlg)                                             \
    d("LArCheatingCosmicRayShowerMatching",     CheatingCosmicRayShowerMatchingAlg)                                             \
    d("LArCheatingNeutrinoCreation",            CheatingNeutrinoCreationAlgorithm)                                              \
    d("LArCheatingNeutrinoDaughterVertices",    CheatingNeutrinoDaughterVerticesAlgorithm)                                      \
    d("LArCheatingPfoCharacterisation",         CheatingPfoCharacterisationAlgorithm)                                           \
    d("LArCheatingPfoCreation",                 CheatingPfoCreationAlgorithm)                                                   \
    d("LArCheatingVertexCreation",              CheatingVertexCreationAlgorithm)                                                \
    d("LArPcaShowerParticleBuilding",           PcaShowerParticleBuildingAlgorithm)                                             \
    d("LArShowerParticleBuilding",              ShowerParticleBuildingAlgorithm)                                                \
    d("LArTrackParticleBuilding",               TrackParticleBuildingAlgorithm)                                                 \
    d("LArStitching",                           StitchingAlgorithm)                                                             \
    d("LArCosmicRayTrackRecovery",              CosmicRayTrackRecoveryAlgorithm)                                                \
    d("LArCosmicRayVertexBuilding",             CosmicRayVertexBuildingAlgorithm)                                               \
    d("LArNeutrinoCreation",                    NeutrinoCreationAlgorithm)                                                      \
    d("LArNeutrinoDaughterVertices",            NeutrinoDaughterVerticesAlgorithm)                                              \
    d("LArNeutrinoHierarchy",                   NeutrinoHierarchyAlgorithm)                                                     \
    d("LArNeutrinoProperties",                  NeutrinoPropertiesAlgorithm)                                                    \
    d("LArDeltaRayIdentification",              DeltaRayIdentificationAlgorithm)                                                \
    d("LArDeltaRayMatching",                    DeltaRayMatchingAlgorithm)                                                      \
    d("LArCosmicRayShowerMatching",             CosmicRayShowerMatchingAlgorithm)                                               \
    d("LArCosmicRayTrackMatching",              CosmicRayTrackMatchingAlgorithm)                                                \
    d("LArThreeDHitCreation",                   ThreeDHitCreationAlgorithm)                                                     \
    d("LArThreeDLongitudinalTracks",            ThreeDLongitudinalTracksAlgorithm)                                              \
    d("LArSlidingConePfoMopUp",                 SlidingConePfoMopUpAlgorithm)                                                   \
    d("LArShowerPfoMopUp",                      ShowerPfoMopUpAlgorithm)                                                        \
    d("LArVertexBasedPfoMopUp",                 VertexBasedPfoMopUpAlgorithm)                                                   \
    d("LArParticleRecovery",                    ParticleRecoveryAlgorithm)                                                      \
    d("LArVertexBasedPfoRecovery",              VertexBasedPfoRecoveryAlgorithm)                                                \
    d("LArThreeDRemnants",                      ThreeDRemnantsAlgorithm)                                                        \
    d("LArThreeDShowers",                       ThreeDShowersAlgorithm)                                                         \
    d("LArThreeDTrackFragments",                ThreeDTrackFragmentsAlgorithm)                                                  \
    d("LArThreeDTransverseTracks",              ThreeDTransverseTracksAlgorithm)                                                \
    d("LArClusterCharacterisation",             ClusterCharacterisationAlgorithm)                                               \
    d("LArPfoCharacterisation",                 PfoCharacterisationAlgorithm)                                                   \
    d("LArShowerGrowing",                       ShowerGrowingAlgorithm)                                                         \
    d("LArCrossGapsAssociation",                CrossGapsAssociationAlgorithm)                                                  \
    d("LArCrossGapsExtension",                  CrossGapsExtensionAlgorithm)                                                    \
    d("LArLongitudinalAssociation",             LongitudinalAssociationAlgorithm)                                               \
    d("LArLongitudinalExtension",               LongitudinalExtensionAlgorithm)                                                 \
    d("LArSimpleClusterGrowing",                SimpleClusterGrowingAlgorithm)                                                  \
    d("LArSimpleClusterMerging",                SimpleClusterMergingAlgorithm)                                                  \
    d("LArTransverseAssociation",               TransverseAssociationAlgorithm)                                                 \
    d("LArTransverseExtension",                 TransverseExtensionAlgorithm)                                                   \
    d("LArSimpleClusterCreation",               SimpleClusterCreationAlgorithm)                                                 \
    d("LArTrackClusterCreation",                TrackClusterCreationAlgorithm)                                                  \
    d("LArClusteringParent",                    ClusteringParentAlgorithm)                                                      \
    d("LArBoundedClusterMopUp",                 BoundedClusterMopUpAlgorithm)                                                   \
    d("LArConeClusterMopUp",                    ConeClusterMopUpAlgorithm)                                                      \
    d("LArIsolatedClusterMopUp",                IsolatedClusterMopUpAlgorithm)                                                  \
    d("LArNearbyClusterMopUp",                  NearbyClusterMopUpAlgorithm)                                                    \
    d("LArSlidingConeClusterMopUp",             SlidingConeClusterMopUpAlgorithm)                                               \
    d("LArCosmicRayExtension",                  CosmicRayExtensionAlgorithm)                                                    \
    d("LArCosmicRaySplitting",                  CosmicRaySplittingAlgorithm)                                                    \
    d("LArDeltaRayExtension",                   DeltaRayExtensionAlgorithm)                                                     \
    d("LArDeltaRayGrowing",                     DeltaRayGrowingAlgorithm)                                                       \
    d("LArBranchSplitting",                     BranchSplittingAlgorithm)                                                       \
    d("LArCrossedTrackSplitting",               CrossedTrackSplittingAlgorithm)                                                 \
    d("LArDeltaRaySplitting",                   DeltaRaySplittingAlgorithm)                                                     \
    d("LArKinkSplitting",                       KinkSplittingAlgorithm)                                                         \
    d("LArLayerSplitting",                      LayerSplittingAlgorithm)                                                        \
    d("LArOvershootSplitting",                  OvershootSplittingAlgorithm)                                                    \
    d("LArTrackConsolidation",                  TrackConsolidationAlgorithm)                                                    \
    d("LArVertexSplitting",                     VertexSplittingAlgorithm)                                                       \
    d("LArTwoDParticleCreation",                TwoDParticleCreationAlgorithm)                                                  \
    d("LArListChanging",                        ListChangingAlgorithm)                                                          \
    d("LArListDeletion",                        ListDeletionAlgorithm)                                                          \
    d("LArListDissolution",                     ListDissolutionAlgorithm)                                                       \
    d("LArListMerging",                         ListMergingAlgorithm)                                                           \
    d("LArListMoving",                          ListMovingAlgorithm)                                                            \
    d("LArListPreparation",                     ListPreparationAlgorithm)                                                       \
    d("LArNeutrinoParent",                      NeutrinoParentAlgorithm)                                                        \
    d("LArCandidateVertexCreation",             CandidateVertexCreationAlgorithm)                                               \
    d("LArEnergyKickVertexSelection",           EnergyKickVertexSelectionAlgorithm)                                             \
    d("LArHitAngleVertexSelection",             HitAngleVertexSelectionAlgorithm)                                               \
    d("LArSVMVertexSelection",                  SVMVertexSelectionAlgorithm)

#define LAR_ALGORITHM_TOOL_LIST(d)                                                                                              \
    d("LArStitchingCosmicRayMerging",           StitchingCosmicRayMergingTool)                                                  \
    d("LArStitchingObjectCreation",             StitchingObjectCreationTool)                                                    \
    d("LArCheatingEventSlicing",                CheatingEventSlicingTool)                                                       \
    d("LArBranchAssociatedPfos",                BranchAssociatedPfosTool)                                                       \
    d("LArEndAssociatedPfos",                   EndAssociatedPfosTool)                                                          \
    d("LArEventSlicing",                        EventSlicingTool)                                                               \
    d("LArVertexAssociatedPfos",                VertexAssociatedPfosTool)                                                       \
    d("LArClearShowers",                        ClearShowersTool)                                                               \
    d("LArShowerTensorVisualization",           ShowerTensorVisualizationTool)                                                  \
    d("LArSimpleShowers",                       SimpleShowersTool)                                                              \
    d("LArSplitShowers",                        SplitShowersTool)                                                               \
    d("LArClearTrackFragments",                 ClearTrackFragmentsTool)                                                        \
    d("LArClearLongitudinalTrackHits",          ClearLongitudinalTrackHitsTool)                                                 \
    d("LArClearTransverseTrackHits",            ClearTransverseTrackHitsTool)                                                   \
    d("LArDeltaRayShowerHits",                  DeltaRayShowerHitsTool)                                                         \
    d("LArMultiValuedLongitudinalTrackHits",    MultiValuedLongitudinalTrackHitsTool)                                           \
    d("LArMultiValuedTransverseTrackHits",      MultiValuedTransverseTrackHitsTool)                                             \
    d("LArThreeViewShowerHits",                 ThreeViewShowerHitsTool)                                                        \
    d("LArTwoViewShowerHits",                   TwoViewShowerHitsTool)                                                          \
    d("LArClearLongitudinalTracks",             ClearLongitudinalTracksTool)                                                    \
    d("LArMatchedEndPoints",                    MatchedEndPointsTool)                                                           \
    d("LArClearRemnants",                       ClearRemnantsTool)                                                              \
    d("LArConnectedRemnants",                   ConnectedRemnantsTool)                                                          \
    d("LArMopUpRemnants",                       MopUpRemnantsTool)                                                              \
    d("LArClearTracks",                         ClearTracksTool)                                                                \
    d("LArLongTracks",                          LongTracksTool)                                                                 \
    d("LArMissingTrack",                        MissingTrackTool)                                                               \
    d("LArMissingTrackSegment",                 MissingTrackSegmentTool)                                                        \
    d("LArOvershootTracks",                     OvershootTracksTool)                                                            \
    d("LArTracksCrossingGaps",                  TracksCrossingGapsTool)                                                         \
    d("LArTrackSplitting",                      TrackSplittingTool)                                                             \
    d("LArTransverseTensorVisualization",       TransverseTensorVisualizationTool)                                              \
    d("LArUndershootTracks",                    UndershootTracksTool)                                                           \
    d("LArEnergyKickFeature",                   EnergyKickFeatureTool)                                                          \
    d("LArGlobalAsymmetryFeature",              GlobalAsymmetryFeatureTool)                                                     \
    d("LArLocalAsymmetryFeature",               LocalAsymmetryFeatureTool)                                                      \
    d("LArRPhiFeature",                         RPhiFeatureTool)                                                                \
    d("LArShowerAsymmetryFeature",              ShowerAsymmetryFeatureTool)

#define LAR_PARTICLE_ID_LIST(d)                                                                                                 \
    d("LArMuonId",                              LArParticleIdPlugins::LArMuonId)

#define FACTORY Factory

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

#define LAR_CONTENT_CREATE_ALGORITHM_FACTORY(a, b)                                                                              \
class b##FACTORY : public pandora::AlgorithmFactory                                                                             \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::Algorithm *CreateAlgorithm() const {return new b;};                                                                \
};

LAR_ALGORITHM_LIST(LAR_CONTENT_CREATE_ALGORITHM_FACTORY)

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY(a, b)                                                                         \
class b##FACTORY : public pandora::AlgorithmToolFactory                                                                         \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::AlgorithmTool *CreateAlgorithmTool() const {return new b;};                                                        \
};

LAR_ALGORITHM_TOOL_LIST(LAR_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY)

} // namespace lar_content

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_REGISTER_ALGORITHM(a, b)                                                                                    \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmFactory(pandora, a, new lar_content::b##FACTORY));        \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

#define LAR_CONTENT_REGISTER_ALGORITHM_TOOL(a, b)                                                                               \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmToolFactory(pandora, a, new lar_content::b##FACTORY));    \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LArContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(LAR_CONTENT_REGISTER_ALGORITHM);
    LAR_ALGORITHM_TOOL_LIST(LAR_CONTENT_REGISTER_ALGORITHM_TOOL);
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_REGISTER_PARTICLE_ID(a, b)                                                                                  \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterParticleIdPlugin(pandora, a, new lar_content::b));                 \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LArContent::RegisterBasicPlugins(const pandora::Pandora &pandora)
{
    LAR_PARTICLE_ID_LIST(LAR_CONTENT_REGISTER_PARTICLE_ID);
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LArContent::SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, lar_content::LArPseudoLayerPlugin *const pLArPseudoLayerPlugin)
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, pLArPseudoLayerPlugin));
    return lar_content::LArGeometryHelper::SetLArPseudoLayerPlugin(pandora, pLArPseudoLayerPlugin);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode LArContent::SetLArTransformationPlugin(const pandora::Pandora &pandora, lar_content::LArTransformationPlugin *const pLArTransformationPlugin)
{
    return lar_content::LArGeometryHelper::SetLArTransformationPlugin(pandora, pLArTransformationPlugin);
}
