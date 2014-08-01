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
#include "LArHelpers/LArPfoHelper.h"

#include "LArMonitoring/EventDisplayAlgorithm.h"
#include "LArMonitoring/NtupleWritingAlgorithm.h"
#include "LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "LArMonitoring/VisualMonitoringAlgorithm.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"
#include "LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"
#include "LArThreeDReco/LArEventBuilding/NeutrinoBuildingAlgorithm.h"
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
#include "LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"
#include "LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"
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

#include "LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.h"
#include "LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"
#include "LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"
#include "LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"
#include "LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"
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
        d("LArCosmicRayIdentification",             lar::CosmicRayIdentificationAlgorithm::Factory)                             \
        d("LArCosmicRayTrackMatching",              lar::CosmicRayTrackMatchingAlgorithm::Factory)                              \
        d("LArNeutrinoBuilding",                    lar::NeutrinoBuildingAlgorithm::Factory)                                    \
        d("LArDeltaRayIdentification",              lar::DeltaRayIdentificationAlgorithm::Factory)                              \
        d("LArDeltaRayMatching",                    lar::DeltaRayMatchingAlgorithm::Factory)                                    \
        d("LArThreeDHitCreation",                   lar::ThreeDHitCreationAlgorithm::Factory)                                   \
        d("LArThreeDLongitudinalTracks",            lar::ThreeDLongitudinalTracksAlgorithm::Factory)                            \
        d("LArThreeDRemnants",                      lar::ThreeDRemnantsAlgorithm::Factory)                                      \
        d("LArThreeDShowers",                       lar::ThreeDShowersAlgorithm::Factory)                                       \
        d("LArThreeDTrackFragments",                lar::ThreeDTrackFragmentsAlgorithm::Factory)                                \
        d("LArThreeDTransverseTracks",              lar::ThreeDTransverseTracksAlgorithm::Factory)                              \
        d("LArLongitudinalAssociation",             lar::LongitudinalAssociationAlgorithm::Factory)                             \
        d("LArLongitudinalExtension",               lar::LongitudinalExtensionAlgorithm::Factory)                               \
        d("LArSimpleClusterMerging",                lar::SimpleClusterMergingAlgorithm::Factory)                                \
        d("LArTransverseAssociation",               lar::TransverseAssociationAlgorithm::Factory)                               \
        d("LArTransverseExtension",                 lar::TransverseExtensionAlgorithm::Factory)                                 \
        d("LArSimpleClusterCreation",               lar::SimpleClusterCreationAlgorithm::Factory)                               \
        d("LArTrackClusterCreation",                lar::TrackClusterCreationAlgorithm::Factory)                                \
        d("LArClusteringParent",                    lar::ClusteringParentAlgorithm::Factory)                                    \
        d("LArBoundedClusterMerging",               lar::BoundedClusterMergingAlgorithm::Factory)                               \
        d("LArConeBasedMerging",                    lar::ConeBasedMergingAlgorithm::Factory)                                    \
        d("LArIsolatedHitMerging",                  lar::IsolatedHitMergingAlgorithm::Factory)                                  \
        d("LArCosmicRayExtension",                  lar::CosmicRayExtensionAlgorithm::Factory)                                  \
        d("LArCosmicRaySplitting",                  lar::CosmicRaySplittingAlgorithm::Factory)                                  \
        d("LArDeltaRayExtension",                   lar::DeltaRayExtensionAlgorithm::Factory)                                   \
        d("LArDeltaRayGrowing",                     lar::DeltaRayGrowingAlgorithm::Factory)                                     \
        d("LArBranchSplitting",                     lar::BranchSplittingAlgorithm::Factory)                                     \
        d("LArCrossedTrackSplitting",               lar::CrossedTrackSplittingAlgorithm::Factory)                               \
        d("LArDeltaRaySplitting",                   lar::DeltaRaySplittingAlgorithm::Factory)                                   \
        d("LArKinkSplitting",                       lar::KinkSplittingAlgorithm::Factory)                                       \
        d("LArLayerSplitting",                      lar::LayerSplittingAlgorithm::Factory)                                      \
        d("LArTrackConsolidation",                  lar::TrackConsolidationAlgorithm::Factory)                                  \
        d("LArVertexSplitting",                     lar::VertexSplittingAlgorithm::Factory)                                     \
        d("LArClusterCharacterisation",             lar::ClusterCharacterisationAlgorithm::Factory)                             \
        d("LArTwoDParticleCreationAlgorithm",       lar::TwoDParticleCreationAlgorithm::Factory)                                \
        d("LArListChanging",                        lar::ListChangingAlgorithm::Factory)                                        \
        d("LArListDissolution",                     lar::ListDissolutionAlgorithm::Factory)                                     \
        d("LArListMerging",                         lar::ListMergingAlgorithm::Factory)                                         \
        d("LArListPreparation",                     lar::ListPreparationAlgorithm::Factory)

    #define LAR_ALGORITHM_TOOL_LIST(d)                                                                                          \
        d("LArClearShowers",                        lar::ClearShowersTool::Factory)                                             \
        d("LArShowerTensorVisualization",           lar::ShowerTensorVisualizationTool::Factory)                                \
        d("LArSimpleShowers",                       lar::SimpleShowersTool::Factory)                                            \
        d("LArSplitShowers",                        lar::SplitShowersTool::Factory)                                             \
        d("LArClearTrackFragments",                 lar::ClearTrackFragmentsTool::Factory)                                      \
        d("LArClearLongitudinalTrackHits",          lar::ClearLongitudinalTrackHitsTool::Factory)                               \
        d("LArClearTransverseTrackHits",            lar::ClearTransverseTrackHitsTool::Factory)                                 \
        d("LArDeltaRayShowerHits",                  lar::DeltaRayShowerHitsTool::Factory)                                       \
        d("LArMultiValuedLongitudinalTrackHits",    lar::MultiValuedLongitudinalTrackHitsTool::Factory)                         \
        d("LArMultiValuedTransverseTrackHits",      lar::MultiValuedTransverseTrackHitsTool::Factory)                           \
        d("LArThreeViewShowerHits",                 lar::ThreeViewShowerHitsTool::Factory)                                      \
        d("LArTwoViewShowerHits",                   lar::TwoViewShowerHitsTool::Factory)                                        \
        d("LArClearLongitudinalTracks",             lar::ClearLongitudinalTracksTool::Factory)                                  \
        d("LArMatchedEndPoints",                    lar::MatchedEndPointsTool::Factory)                                         \
        d("LArClearRemnants",                       lar::ClearRemnantsTool::Factory)                                            \
        d("LArClearTracks",                         lar::ClearTracksTool::Factory)                                              \
        d("LArLongTracks",                          lar::LongTracksTool::Factory)                                               \
        d("LArMissingTrack",                        lar::MissingTrackTool::Factory)                                             \
        d("LArMissingTrackSegment",                 lar::MissingTrackSegmentTool::Factory)                                      \
        d("LArOvershootTracks",                     lar::OvershootTracksTool::Factory)                                          \
        d("LArTrackSplitting",                      lar::TrackSplittingTool::Factory)                                           \
        d("LArTransverseTensorVisualization",       lar::TransverseTensorVisualizationTool::Factory)                            \
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
        d("LArPfoHelper",                           &lar::LArPfoHelper::ReadSettings)

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

inline pandora::StatusCode LArContent::RegisterResetFunctions(pandora::Pandora &/*pandora*/)
{
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
