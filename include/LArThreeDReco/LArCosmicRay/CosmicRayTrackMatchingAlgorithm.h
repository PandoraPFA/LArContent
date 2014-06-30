/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray track matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  CosmicRayTrackMatchingAlgorithm class
 */
class CosmicRayTrackMatchingAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pClusterU the cluster in the U view
         *  @param  pClusterV the cluster in the V view
         *  @param  pClusterW the cluster in the W view
         */
        Particle(const pandora::Cluster *pClusterU, const pandora::Cluster *pClusterV, const pandora::Cluster *pClusterW);

        const pandora::Cluster  *m_pClusterU;    ///< Address of cluster in U view
        const pandora::Cluster  *m_pClusterV;    ///< Address of cluster in V view
        const pandora::Cluster  *m_pClusterW;    ///< Address of cluster in W view
    };

    typedef std::vector<Particle> ParticleList;

    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterAssociationMap;

    /**
     *  @brief Get a vector of available clusters
     *
     *  @param inputClusterListName the input name of the cluster list
     *  @param clusterVector the output vector of available clusters
     */
    pandora::StatusCode GetAvailableClusters(const std::string inputClusterListName, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief Select a set of clusters judged to be clean
     *
     *  @param inputVector the input vector of all available clusters
     *  @param outputVector the output vector of clean clusters
     */
    void SelectCleanClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;

    /**
     *  @brief Match sets of clusters from two views
     *
     *  @param clusterVector1 the vector of clusters from the first view
     *  @param clusterVector2 the vector of clusters from the second view
     *  @param matchedClusters12 the map of cluster matches
     */
    void MatchClusters(const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2,
        CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12) const;

    /**
     *  @brief Match a pair of clusters from two views
     *
     *  @param pCluster1 the first cluster
     *  @param pCluster2 the second cluster
     *
     *  @return boolean
     */
    bool MatchClusters(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief Match clusters from three views and form into particles
     *
     *  @param matchedClusters12 the map of matches between the view 1 and view 2
     *  @param matchedClusters23 the map of matches between the view 2 and view 3
     *  @param matchedClusters31 the map of matches between the view 3 and view 1
     *  @param particleList the output list of particles
     */
    void MatchThreeViews(const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters23,
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters31,
        ParticleList &particleList) const;

    /**
     *  @brief Match clusters from two views and form into particles
     *
     *  @param matchedClusters12 the map of matches between the view 1 and view 2
     *  @param matchedClusters23 the map of matches between the view 2 and view 3
     *  @param matchedClusters31 the map of matches between the view 3 and view 1
     *  @param particleList the output list of particles
     */
    void MatchTwoViews(const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters23,
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters31,
        ParticleList &particleList) const;

    /**
     *  @brief Match clusters from two views and form into particles
     *
     *  @param matchedClusters12 the map of matches between the first and second views
     *  @param particleList the output list of particles
     */
    void MatchTwoViews(const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
        ParticleList &particleList) const;

    /**
     *  @brief Remove ambiguities between candidate particles
     *
     *  @param inputList the input list of particles
     *  @param outputList the output list of particles
     */
    void ResolveAmbiguities(const ParticleList &inputList, ParticleList &outputList) const;

    /**
     *  @brief Build PFO objects from candidate particles
     *
     *  @param particleList the input list of particles
     */
    void BuildParticles(const ParticleList &particleList);

    /**
     *  @brief Check that three clusters have a consistent 3D position
     *
     *  @param pCluster1 the cluster from the first view
     *  @param pCluster2 the cluster from the second view
     *  @param pCluster3 the cluster from the third view
     *
     *  @return boolean
     */
    bool CheckMatchedClusters3D(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::Cluster *const pCluster3) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string    m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string    m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string    m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string    m_outputPfoListName;            ///< The name of the output PFO list

    float          m_clusterMinLength;             ///< minimum length of clusters for this algorithm
    float          m_vtxXOverlap;                  ///< requirement on X overlap of start/end positions
    float          m_minXOverlap;                  ///< requirement on minimum X overlap for associated clusters
    float          m_minXOverlapFraction;          ///< requirement on minimum X overlap fraction for associated clusters
    float          m_maxDisplacement;              ///< requirement on 3D consistency checks
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
