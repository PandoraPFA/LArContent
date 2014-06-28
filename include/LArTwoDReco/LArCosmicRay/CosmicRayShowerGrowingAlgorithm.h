/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray shower growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SHOWER_GROWING_ALGORITHM_H
#define LAR_COSMIC_RAY_SHOWER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRayShowerGrowingAlgorithm class
 */
class CosmicRayShowerGrowingAlgorithm : public CosmicRayGrowingAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    void GetListOfSeedClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Use primary Pfos to select seed clusters
     *
     *  @param  pfoVector the input vector of primary Pfos
     *  @param  hitType the cluster hit type
     *  @param  seedClusters the output vector of seed clusters
     */
    void SelectPrimaryPfoSeeds(const pandora::PfoVector &pfoVector, const pandora::HitType hitType,
        pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Use secondary Pfos to select seed clusters
     *
     *  @param  pfoVector the input vector of secondary Pfos
     *  @param  hitType the cluster hit type
     *  @param  seedClusters the output vector of seed clusters
     */
    void SelectSecondaryPfoSeeds(const pandora::PfoVector &pfoVector, const pandora::HitType hitType,
        pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Use available clusters to select seed clusters
     *
     *  @param  clusterVector the input vector of clean clusters
     *  @param  pfoVector the input vector of primary Pfos
     *  @param  hitType the cluster hit type
     *  @param  seedClusters the output vector of seed clusters
     */
    void SelectClusterSeeds(const pandora::ClusterVector &clusterVector, const pandora::PfoVector &pfoVector,
        const pandora::HitType hitType, pandora::ClusterVector &seedClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool         m_growPfos;                     ///< Grow pfo clusters
    bool         m_growClusters;                 ///< Grow available clusters
    float        m_maxSeedClusterLength;         ///< The maximum length of primary Pfo's that can be taken as seed clusters
    float        m_maxSeedClusterDisplacement;   ///< The maximum displacement of a possible seed cluster from a primary Pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerGrowingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_GROWING_ALGORITHM_H
