/**
 *  @file   LArContent/include/LArContent.h
 * 
 *  @brief  Header file detailing content for use with liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "ClusterAssociation/ClusterAssociationAlgorithm.h"
#include "ClusterAssociation/ClusterExtensionAlgorithm.h"
#include "ClusterAssociation/IsolatedHitMergingAlgorithm.h"
#include "Clustering/ClusterCreationAlgorithm.h"
#include "Clustering/ClusteringParentAlgorithm.h"
#include "Clustering/RemnantClusteringAlgorithm.h"
#include "Clustering/TransverseClusteringAlgorithm.h"
#include "ClusterSeedAssociation/BoundedClusterMergingAlgorithm.h"
#include "ClusterSeedAssociation/ConeBasedMergingAlgorithm.h"
#include "ClusterSeedAssociation/ParallelClusterMergingAlgorithm.h"
#include "ClusterSplitting/KinkSplittingAlgorithm.h"
#include "Monitoring/EventDisplayAlgorithm.h"
#include "Monitoring/NtupleWritingAlgorithm.h"
#include "Monitoring/VisualMonitoringAlgorithm.h"
#include "Reclustering/ShowerMipSeparationAlgorithm.h"
#include "ThreeDSeed/ThreeDTrackSegmentsAlgorithm.h"
#include "TwoDSeed/SeedBranchGrowingAlgorithm.h"
#include "TwoDSeed/SeedConsolidationAlgorithm.h"
#include "TwoDSeed/SeedFindingAlgorithm.h"
#include "TwoDSeed/SeedLengthGrowingAlgorithm.h"
#include "TwoDSeed/SeedRelegationAlgorithm.h"
#include "TwoDSeed/VertexSeedFindingAlgorithm.h"
#include "Utility/EventPreparationAlgorithm.h"
#include "Utility/TwoDPreparationAlgorithm.h"
#include "Vertex/VertexFindingAlgorithm.h"
#include "Vertex/VertexSplittingAlgorithm.h"

#include "Helpers/LArClusterHelper.h"
#include "Helpers/LArGeometryHelper.h"
#include "Helpers/LArParticleIdHelper.h"
#include "Helpers/LArPointingClusterHelper.h"
#include "Helpers/LArVertexHelper.h"

#include "LArBFieldCalculator.h"
#include "LArPseudoLayerCalculator.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArClusterAssociation",                  lar::ClusterAssociationAlgorithm::Factory)                                  \
        d("LArClusterExtension",                    lar::ClusterExtensionAlgorithm::Factory)                                    \
        d("LArIsolatedHitMerging",                  lar::IsolatedHitMergingAlgorithm::Factory)                                  \
        d("LArClusterCreation",                     lar::ClusterCreationAlgorithm::Factory)                                     \
        d("LArClusteringParent",                    lar::ClusteringParentAlgorithm::Factory)                                    \
        d("LArRemnantClustering",                   lar::RemnantClusteringAlgorithm::Factory)                                   \
        d("LArTransverseClustering",                lar::TransverseClusteringAlgorithm::Factory)                                \
        d("LArBoundedClusterMerging",               lar::BoundedClusterMergingAlgorithm::Factory)                               \
        d("LArConeBasedMerging",                    lar::ConeBasedMergingAlgorithm::Factory)                                    \
        d("LArParallelClusterMerging",              lar::ParallelClusterMergingAlgorithm::Factory)                              \
        d("LArKinkSplitting",                       lar::KinkSplittingAlgorithm::Factory)                                       \
        d("LArEventDisplay",                        lar::EventDisplayAlgorithm::Factory)                                        \
        d("LArNtupleWriting",                       lar::NtupleWritingAlgorithm::Factory)                                       \
        d("LArVisualMonitoring",                    lar::VisualMonitoringAlgorithm::Factory)                                    \
        d("LArShowerMipSeparation",                 lar::ShowerMipSeparationAlgorithm::Factory)                                 \
        d("LArThreeDTrackSegments",                 lar::ThreeDTrackSegmentsAlgorithm::Factory)                                 \
        d("LArSeedBranchGrowing",                   lar::SeedBranchGrowingAlgorithm::Factory)                                   \
        d("LArSeedConsolidation",                   lar::SeedConsolidationAlgorithm::Factory)                                   \
        d("LArSeedFinding",                         lar::SeedFindingAlgorithm::Factory)                                         \
        d("LArSeedLengthGrowing",                   lar::SeedLengthGrowingAlgorithm::Factory)                                   \
        d("LArSeedRelegation",                      lar::SeedRelegationAlgorithm::Factory)                                      \
        d("LArVertexSeedFinding",                   lar::VertexSeedFindingAlgorithm::Factory)                                   \
        d("LArEventPreparation",                    lar::EventPreparationAlgorithm::Factory)                                    \
        d("LArTwoDPreparation",                     lar::TwoDPreparationAlgorithm::Factory)                                     \
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
