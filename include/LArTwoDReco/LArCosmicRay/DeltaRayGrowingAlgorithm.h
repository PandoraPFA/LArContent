/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the delta ray growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_GROWING_ALGORITHM_H
#define LAR_DELTA_RAY_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  DeltaRayGrowingAlgorithm class
 */
class DeltaRayGrowingAlgorithm : public CosmicRayGrowingAlgorithm
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
    
    float        m_maxPrimaryClusterDisplacement;   ///< The maximum displacement of a possible seed cluster from a primary Pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayGrowingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_DELTA_RAY_GROWING_ALGORITHM_H
