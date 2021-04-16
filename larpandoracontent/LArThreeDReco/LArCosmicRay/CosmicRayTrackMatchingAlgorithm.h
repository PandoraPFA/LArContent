/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray track matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayBaseMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayTrackMatchingAlgorithm class
 */
class CosmicRayTrackMatchingAlgorithm : public CosmicRayBaseMatchingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRayTrackMatchingAlgorithm();

private:
    void SelectCleanClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;
    bool MatchClusters(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;
    bool CheckMatchedClusters3D(
        const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3) const;
    void SetPfoParameters(const Particle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_clusterMinLength;    ///< minimum length of clusters for this algorithm
    float m_vtxXOverlap;         ///< requirement on X overlap of start/end positions
    float m_minXOverlap;         ///< requirement on minimum X overlap for associated clusters
    float m_minXOverlapFraction; ///< requirement on minimum X overlap fraction for associated clusters
    float m_maxDisplacement;     ///< requirement on 3D consistency checks
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
