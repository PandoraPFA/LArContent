/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayGrowingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_GROWING_ALGORITHM_H
#define LAR_COSMIC_RAY_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRayGrowingAlgorithm class
 */
class CosmicRayGrowingAlgorithm : public ClusterGrowingAlgorithm
{
protected:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &cleanClusters) const;
    virtual void GetListOfSeedClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &seedClusters) const = 0;

    /**
     *  @brief  Get a vector of Pfos from an input Pfo list name
     *
     *  @param  inputPfoListName the input Pfo list name
     *  @param  pfoVector the output vector of Pfos
     */
    void GetPfos(const std::string inputPfoListName, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Get a list of clusters of a particular hit type from a vector of pfos
     *
     *  @param  pfos the input vector of Pfos
     *  @param  hitType the cluster hit type
     *  @param  clusters the output list of clusters
     */
    void GetPfoClusters(const pandora::PfoVector &pfos, const pandora::HitType hitType, pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Get a list of clusters of a particular hit type from a given pfo
     *
     *  @param  pfo the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  clusters the output list of clusters
     */
    void GetPfoClusters(const pandora::ParticleFlowObject *pfo, const pandora::HitType hitType, pandora::ClusterList &clusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_primaryPfoListName;            ///< The primary Pfo list name
    std::string m_secondaryPfoListName;          ///< The secondary Pfo list name

    unsigned int m_minCaloHitsPerCluster;        ///< The min number of calo hits per candidate cluster
};

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_GROWING_ALGORITHM_H
