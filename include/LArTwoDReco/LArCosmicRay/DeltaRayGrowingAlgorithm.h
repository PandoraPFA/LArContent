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

#include "LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  DeltaRayGrowingAlgorithm class
 */
class DeltaRayGrowingAlgorithm : public ClusterGrowingAlgorithm
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
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &cleanClusters) const;
    void GetListOfSeedClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Get a vector of Pfos from an input Pfo list name
     *
     *  @param  inputPfoListName the input Pfo list name
     *  @param  pfoVector the output vector of Pfos
     */
    void GetPfos(const std::string inputPfoListName, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Use available clusters to select seed clusters
     *
     *  @param  clusterVector the input vector of clean clusters
     *  @param  pfoVector the input vector of primary Pfos
     *  @param  hitType the cluster hit type
     *  @param  seedClusters the output vector of seed clusters
     */
    void SelectSeedClusters(const pandora::ClusterVector &clusterVector, const pandora::PfoVector &pfoVector,
        const pandora::HitType hitType, pandora::ClusterVector &seedClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputPfoListName;                 ///< The input Pfo list name

    unsigned int m_minCaloHitsPerCluster;           ///< The minimum number of calo hits per candidate cluster
    unsigned int m_minSeedClusterCaloHits;          ///< The minimum number of calo hits for seed clusters
    float        m_maxSeedClusterDisplacement;      ///< The maximum displacement of a seed cluster from a primary track
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayGrowingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_DELTA_RAY_GROWING_ALGORITHM_H
