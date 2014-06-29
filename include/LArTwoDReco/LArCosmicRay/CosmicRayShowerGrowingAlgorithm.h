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

#include "LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRayShowerGrowingAlgorithm class
 */
class CosmicRayShowerGrowingAlgorithm : public ClusterGrowingAlgorithm
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
     *  @brief  Select seed cluster for growing
     *
     *  @param  clusterVector the input vector of clusters
     *  @param  pfoVector the input vector of pPfos
     *  @param  hitType the cluster hit type
     *  @param  seedClusters the output vector of seed clusters
     */
    void SelectSeedClusters(const pandora::ClusterVector &inputClusters, const pandora::PfoVector &pfoVector, 
        const pandora::HitType &hitType, pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Determine whether a Pfo is a possible delta ray
     *
     *  @param  pPfo the input particle flow object
     *  @param  pfoVector the input vector of particle flow vectors
     *
     *  @return boolean
     */
    bool IsPossibleDeltaRay(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoVector &pfoVector) const;


    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_inputPfoListName;                ///< The primary Pfo list name

    unsigned int    m_minCaloHitsPerCluster;           ///< The overall minimum number of calo hits for cluster
    unsigned int    m_minSeedClusterCaloHits;          ///< The minimum number of calo hits in seed clusters
    float           m_maxSeedClusterLength;            ///< The maximum length of a unassociated shower clusters
    float           m_maxSeedClusterDisplacement;      ///< The maximum displacement for associating shower clusters with tracks
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerGrowingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_GROWING_ALGORITHM_H
