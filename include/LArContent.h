/**
 *  @file   LArContent.h
 * 
 *  @brief  Header file detailing content for use with liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "BoundedClusterMergingAlgorithm.h"
#include "ClusterAssociationAlgorithm.h"
#include "ClusterCreationAlgorithm.h"
#include "ClusterExtensionAlgorithm.h"
#include "ClusteringParentAlgorithm.h"
#include "ConeBasedMergingAlgorithm.h"
#include "EventDisplayAlgorithm.h"
#include "EventPreparationAlgorithm.h"
#include "IsolatedHitMergingAlgorithm.h"
#include "KinkSplittingAlgorithm.h"
#include "NtupleWritingAlgorithm.h"
#include "ParallelClusterMergingAlgorithm.h"
#include "RemnantClusteringAlgorithm.h"
#include "SeedBranchGrowingAlgorithm.h"
#include "SeedConsolidationAlgorithm.h"
#include "SeedFindingAlgorithm.h"
#include "SeedLengthGrowingAlgorithm.h"
#include "SeedRelegationAlgorithm.h"
#include "ShowerMipSeparationAlgorithm.h"
#include "ThreeDParticleCreationAlgorithm.h"
#include "ThreeDParticleMatchingAlgorithm.h"
#include "TransverseClusteringAlgorithm.h"
#include "TwoDParticleCreationAlgorithm.h"
#include "TwoDPreparationAlgorithm.h"
#include "VertexFindingAlgorithm.h"
#include "VertexSeedFindingAlgorithm.h"
#include "VisualMonitoringAlgorithm.h"

#include "LArClusterHelper.h"
#include "LArGeometryHelper.h"
#include "LArParticleId.h"
#include "LArPointingClusterHelper.h"
#include "LArVertexHelper.h"

#include "LArBFieldCalculator.h"
#include "LArPseudoLayerCalculator.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArBoundedClusterMerging",               lar::BoundedClusterMergingAlgorithm::Factory)                               \
        d("LArClusterAssociation",                  lar::ClusterAssociationAlgorithm::Factory)                                  \
        d("LArClusterCreation",                     lar::ClusterCreationAlgorithm::Factory)                                     \
        d("LArClusterExtension",                    lar::ClusterExtensionAlgorithm::Factory)                                    \
        d("LArClusteringParent",                    lar::ClusteringParentAlgorithm::Factory)                                    \
        d("LArConeBasedMerging",                    lar::ConeBasedMergingAlgorithm::Factory)                                    \
        d("LArEventDisplay",                        lar::EventDisplayAlgorithm::Factory)                                        \
        d("LArEventPreparation",                    lar::EventPreparationAlgorithm::Factory)                                    \
        d("LArIsolatedHitMerging",                  lar::IsolatedHitMergingAlgorithm::Factory)                                  \
        d("LArKinkSplitting",                       lar::KinkSplittingAlgorithm::Factory)                                       \
        d("LArNtupleWriting",                       lar::NtupleWritingAlgorithm::Factory)                                       \
        d("LArParallelClusterMerging",              lar::ParallelClusterMergingAlgorithm::Factory)                              \
        d("LArRemnantClustering",                   lar::RemnantClusteringAlgorithm::Factory)                                   \
        d("LArSeedBranchGrowing",                   lar::SeedBranchGrowingAlgorithm::Factory)                                   \
        d("LArSeedConsolidation",                   lar::SeedConsolidationAlgorithm::Factory)                                   \
        d("LArSeedFinding",                         lar::SeedFindingAlgorithm::Factory)                                         \
        d("LArSeedLengthGrowing",                   lar::SeedLengthGrowingAlgorithm::Factory)                                   \
        d("LArSeedRelegation",                      lar::SeedRelegationAlgorithm::Factory)                                      \
        d("LArShowerMipSeparation",                 lar::ShowerMipSeparationAlgorithm::Factory)                                 \
        d("LArThreeDParticleCreation",              lar::ThreeDParticleCreationAlgorithm::Factory)                              \
        d("LArThreeDParticleMatching",              lar::ThreeDParticleMatchingAlgorithm::Factory)                              \
        d("LArTransverseClustering",                lar::TransverseClusteringAlgorithm::Factory)                                \
        d("LArTwoDParticleCreation",                lar::TwoDParticleCreationAlgorithm::Factory)                                \
        d("LArTwoDPreparation",                     lar::TwoDPreparationAlgorithm::Factory)                                     \
        d("LArVertexFinding",                       lar::VertexFindingAlgorithm::Factory)                                       \
        d("LArVertexSeedFinding",                   lar::VertexSeedFindingAlgorithm::Factory)                                   \
        d("LArVisualMonitoring",                    lar::VisualMonitoringAlgorithm::Factory)

    #define LAR_PARTICLE_ID_LIST(d)                                                                                             \
        d("LArEmShowerId",                          &lar::LArParticleId::LArEmShowerId)                                         \
        d("LArPhotonId",                            &lar::LArParticleId::LArPhotonId)                                           \
        d("LArElectronId",                          &lar::LArParticleId::LArElectronId)                                         \
        d("LArMuonId",                              &lar::LArParticleId::LArMuonId)

    #define LAR_SETTINGS_LIST(d)                                                                                                \
        d("LArClusterHelper",                       &lar::LArClusterHelper::ReadSettings)                                       \
        d("LArGeometryHelper",                      &lar::LArGeometryHelper::ReadSettings)                                      \
        d("LArParticleId",                          &lar::LArParticleId::ReadSettings)                                          \
        d("LArPointingClusterHelper",               &lar::LArPointingClusterHelper::ReadSettings)                               \
        d("LArVertexHelper",                        &lar::LArVertexHelper::ReadSettings)

    /**
     *  @brief  Register all the fine granularity algorithms with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(pandora::Pandora &pandora);

    /**
     *  @brief  Register all the fine granularity helper functions with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterHelperFunctions(pandora::Pandora &pandora);
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

#endif // #ifndef LAR_CONTENT_H
