/**
 *  @file   LArContent/include/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h
 * 
 *  @brief  Header file for the delta ray matching algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  DeltaRayMatchingAlgorithm class
 */
class DeltaRayMatchingAlgorithm : public pandora::Algorithm
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Cosmic ray to shower matching using shower information from all three views
     * 
     *  @param  clustersU the list of clusters in the u view
     *  @param  clustersV the list of clusters in the v view
     *  @param  clustersW the list of clusters in the w view
     */
    void ThreeViewMatching(const pandora::ClusterVector &clustersU, const pandora::ClusterVector &clustersV, const pandora::ClusterVector &clustersW) const;

    /**
     *  @brief  Cosmic ray to shower matching using shower information from two views
     * 
     *  @param  clusters1 the list of clusters in the first view
     *  @param  clusters2 the list of clusters in the second view
     */
    void TwoViewMatching(const pandora::ClusterVector &clusters1, const pandora::ClusterVector &clusters2) const;

    /**
     *  @brief  Cosmic ray to shower matching using shower information from a single view
     * 
     *  @param  clustersU the list of clusters in the provided view
     */
    void OneViewMatching(const pandora::ClusterVector &clusters) const;

    /**
     *  @brief  Get a vector containing all available input clusters in the provided cluster lists, storing sliding
     *          linear fits in the algorithm fit result cache.
     *
     *  @param  clusterListNames the vector of cluster list names
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetInputClusters(const std::string &clusterListName, pandora::ClusterVector &clusterVector);

    /**
     *  @brief  Look at consistency of a UVW combination
     *
     *  @param  pointer to U view cluster
     *  @param  pointer to V view cluster
     *  @param  pointer to W view cluster
     * 
     *  @return pseudo chi2 of consistency
     */
    float CompareClusterTriplet(pandora::Cluster *const pClusterU, pandora::Cluster *const pClusterV, pandora::Cluster *const pClusterW) const;

    /**
     *  @brief  Return indicative coordinates of cluster at a position x using range indicated
     *
     *  @param  pointer to cluster
     *  @param  x (time) value at which coordinate is returned
     *  @param  minimum x value used to evaluate coordinate
     *  @param  maximum x value used to evaluate coordinate
     */
    float GetCoordinateAtX(pandora::Cluster *const pCluster, const float x, const float xmin, const float xmax) const;

    /**
     *  @brief  Find best PFO to associate a UVW triplet
     *
     *  @param  pointer to U view cluster
     *  @param  pointer to V view cluster
     *  @param  pointer to W view cluster
     *  @param  pointer to best PFO
     *  @param  distance measure to best PFO
     */
    void FindBestCosmicPFO(pandora::Cluster *const pClusterU, pandora::Cluster *const pClusterV, pandora::Cluster *const pClusterW,
        pandora::ParticleFlowObject *&pBestPFO, float &distanceToBestPFO) const;

    /**
     *  @brief  Check that matched cluster is a good sub-cluster of the PFO cluster
     *
     *  @param  pointer to PFO cluster
     *  @param  pointer to matched cluster
     *
     *  @return boolean
     */
    bool IsSubCluster(const pandora::Cluster *const pPFOCluster, const pandora::Cluster *const pSubCluster) const;

    /**
     *  @brief  Create a new pfo using a provided list of clusters and set it to be the daughter of a provided parent pfo
     * 
     *  @param  clusterList the list of clusters
     *  @param  pParentPfo address of the parent pfo
     */
    void CreateDaughterPfo(const pandora::ClusterList &clusterList, pandora::ParticleFlowObject *const pParentPfo) const;

    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingFitResult &GetCachedSlidingFitResult(pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Add a new sliding fit result, for the specified cluster, to the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(pandora::Cluster *const pCluster);

    std::string             m_inputPfoListName;           ///< The input pfo list name
    std::string             m_outputPfoListName;          ///< The output pfo list name for new daughter particles

    std::string             m_inputClusterListNameU;     ///< The input cluster list name for the u view
    std::string             m_inputClusterListNameV;     ///< The input cluster list name for the v view
    std::string             m_inputClusterListNameW;     ///< The input cluster list name for the w view

    float                   m_distanceFor1ViewMatching;   ///< Distance cut for single view matching
    float                   m_distanceFor2ViewMatching;   ///< Distance cut for two view matching
    float                   m_distanceFor3ViewMatching;   ///< Distance cut for three view matching
    float                   m_chi2For3ViewMatching;       ///< Pseudo chi2 cut for three view matching

    unsigned int            m_minCaloHitsPerCluster;      ///< The min number of calo hits per candidate cluster

    unsigned int            m_slidingFitWindow;           ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap m_slidingFitResultMap;        ///< The sliding fit result map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_DELTA_RAY_MATCHING_ALGORITHM_H
