/**
 *  @file   LArContent/include/LArContent.h
 * 
 *  @brief  Header file detailing content for use with liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "LArClusterAssociation/IsolatedHitMergingAlgorithm.h"
#include "LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "LArClustering/ClusterCreationAlgorithm.h"
#include "LArClustering/ClusteringParentAlgorithm.h"
#include "LArClusterSeedAssociation/BoundedClusterMergingAlgorithm.h"
#include "LArClusterSeedAssociation/ConeBasedMergingAlgorithm.h"
#include "LArClusterSeedAssociation/ParallelClusterMergingAlgorithm.h"
#include "LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "LArMonitoring/EventDisplayAlgorithm.h"
#include "LArMonitoring/NtupleWritingAlgorithm.h"
#include "LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "LArMonitoring/VisualMonitoringAlgorithm.h"
#include "LArReclustering/ShowerRebuildingAlgorithm.h"
#include "LArReclustering/TrackSplittingAlgorithm.h"
#include "LArThreeDSeed/ThreeDShowersAlgorithm.h"
#include "LArThreeDSeed/ThreeDLongitudinalTracksAlgorithm.h"
#include "LArThreeDSeed/ThreeDTransverseTracksAlgorithm.h"
#include "LArThreeDSeed/ThreeDPairedTracksAlgorithm.h"
#include "LArTwoDSeed/LengthSeedFindingAlgorithm.h"
#include "LArTwoDSeed/SeedBranchGrowingAlgorithm.h"
#include "LArTwoDSeed/SeedCharacterisationAlgorithm.h"
#include "LArTwoDSeed/SeedConsolidationAlgorithm.h"
#include "LArTwoDSeed/SeedLengthGrowingAlgorithm.h"
#include "LArTwoDSeed/SeedMergingAlgorithm.h"
#include "LArTwoDSeed/SeedRelegationAlgorithm.h"
#include "LArTwoDSeed/VertexSeedFindingAlgorithm.h"
#include "LArTwoDSeed/TwoDParticleCreationAlgorithm.h"
#include "LArUtility/EventPreparationAlgorithm.h"
#include "LArUtility/HelloWorldAlgorithm.h"
#include "LArUtility/ListMergingAlgorithm.h"
#include "LArUtility/TwoDPreparationAlgorithm.h"
#include "LArUtility/ThreeDPreparationAlgorithm.h"
#include "LArVertex/VertexFindingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArCalculators/LArBFieldCalculator.h"
#include "LArCalculators/LArTransformationCalculator.h"
#include "LArCalculators/LArPseudoLayerCalculator.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArIsolatedHitMerging",                  lar::IsolatedHitMergingAlgorithm::Factory)                                  \
        d("LArClusterCreation",                     lar::ClusterCreationAlgorithm::Factory)                                     \
        d("LArClusteringParent",                    lar::ClusteringParentAlgorithm::Factory)                                    \
        d("LArBoundedClusterMerging",               lar::BoundedClusterMergingAlgorithm::Factory)                               \
        d("LArConeBasedMerging",                    lar::ConeBasedMergingAlgorithm::Factory)                                    \
        d("LArParallelClusterMerging",              lar::ParallelClusterMergingAlgorithm::Factory)                              \
        d("LArKinkSplitting",                       lar::KinkSplittingAlgorithm::Factory)                                       \
        d("LArLongitudinalAssociation",             lar::LongitudinalAssociationAlgorithm::Factory)                             \
        d("LArLongitudinalExtension",               lar::LongitudinalExtensionAlgorithm::Factory)                               \
        d("LArEventDisplay",                        lar::EventDisplayAlgorithm::Factory)                                        \
        d("LArNtupleWriting",                       lar::NtupleWritingAlgorithm::Factory)                                       \
        d("LArParticleMonitoring",                  lar::ParticleMonitoringAlgorithm::Factory)                                  \
        d("LArVisualMonitoring",                    lar::VisualMonitoringAlgorithm::Factory)                                    \
        d("LArShowerRebuilding",                    lar::ShowerRebuildingAlgorithm::Factory)                                    \
        d("LArTrackSplitting",                      lar::TrackSplittingAlgorithm::Factory)                                      \
        d("LArThreeDShowers",                       lar::ThreeDShowersAlgorithm::Factory)                                       \
        d("LArThreeDLongitudinalTracks",            lar::ThreeDLongitudinalTracksAlgorithm::Factory)                            \
        d("LArThreeDTransverseTracks",              lar::ThreeDTransverseTracksAlgorithm::Factory)                              \
        d("LArThreeDPairedTracks",                  lar::ThreeDPairedTracksAlgorithm::Factory)                                  \
        d("LArLengthSeedFinding",                   lar::LengthSeedFindingAlgorithm::Factory)                                   \
        d("LArSeedBranchGrowing",                   lar::SeedBranchGrowingAlgorithm::Factory)                                   \
        d("LArSeedCharacterisation",                lar::SeedCharacterisationAlgorithm::Factory)                                \
        d("LArSeedConsolidation",                   lar::SeedConsolidationAlgorithm::Factory)                                   \
        d("LArSeedLengthGrowing",                   lar::SeedLengthGrowingAlgorithm::Factory)                                   \
        d("LArSeedMerging",                         lar::SeedMergingAlgorithm::Factory)                                         \
        d("LArSeedRelegation",                      lar::SeedRelegationAlgorithm::Factory)                                      \
        d("LArVertexSeedFinding",                   lar::VertexSeedFindingAlgorithm::Factory)                                   \
        d("LArTwoDParticleCreationAlgorithm",       lar::TwoDParticleCreationAlgorithm::Factory)                                \
        d("LArTransverseAssociation",               lar::TransverseAssociationAlgorithm::Factory)                               \
        d("LArEventPreparation",                    lar::EventPreparationAlgorithm::Factory)                                    \
        d("LArHelloWorld",                          lar::HelloWorldAlgorithm::Factory)                                          \
        d("LArListMerging",                         lar::ListMergingAlgorithm::Factory)                                         \
        d("LArTwoDPreparation",                     lar::TwoDPreparationAlgorithm::Factory)                                     \
        d("LArThreeDPreparation",                   lar::ThreeDPreparationAlgorithm::Factory)                                   \
        d("LArVertexFinding",                       lar::VertexFindingAlgorithm::Factory)                                       \
        d("LArVertexSplitting",                     lar::VertexSplittingAlgorithm::Factory)

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
        d("LArThreeDHelper",                        &lar::LArThreeDHelper::ReadSettings)                                        \
        d("LArVertexHelper",                        &lar::LArVertexHelper::ReadSettings)

    /**
     *  @brief  Register all the lar content algorithms with pandora
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
     *  @param  pLArBFieldCalculator the address of the lar b field calculator
     */
    static pandora::StatusCode SetLArBFieldCalculator(pandora::Pandora &pandora, lar::LArBFieldCalculator *pLArBFieldCalculator);

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
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterResetFunction(pandora, &lar::LArVertexHelper::Reset));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterResetFunction(pandora, &lar::LArThreeDHelper::Reset));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArBFieldCalculator(pandora::Pandora &pandora, lar::LArBFieldCalculator *pLArBFieldCalculator)
{
    return PandoraApi::SetBFieldCalculator(pandora, pLArBFieldCalculator);
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
