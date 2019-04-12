/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayBaseMatchingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray base matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_BASE_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_BASE_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CosmicRayBaseMatchingAlgorithm class
 */
class CosmicRayBaseMatchingAlgorithm : public pandora::Algorithm
{
protected:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

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
        Particle(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

        const pandora::Cluster  *m_pClusterU;    ///< Address of cluster in U view
        const pandora::Cluster  *m_pClusterV;    ///< Address of cluster in V view
        const pandora::Cluster  *m_pClusterW;    ///< Address of cluster in W view
    };

    typedef std::vector<Particle> ParticleList;
    typedef std::unordered_map<const pandora::Cluster*, pandora::ClusterList> ClusterAssociationMap;

    /**
     *  @brief Select a set of clusters judged to be clean
     *
     *  @param inputVector the input vector of all available clusters
     *  @param outputVector the output vector of clean clusters
     */
    virtual void SelectCleanClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const = 0;

    /**
     *  @brief Match a pair of clusters from two views
     *
     *  @param pCluster1 the first cluster
     *  @param pCluster2 the second cluster
     *
     *  @return boolean
     */
    virtual bool MatchClusters(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const = 0;

    /**
     *  @brief Check that three clusters have a consistent 3D position
     *
     *  @param pCluster1 the cluster from the first view
     *  @param pCluster2 the cluster from the second view
     *  @param pCluster3 the cluster from the third view
     *
     *  @return boolean
     */
    virtual bool CheckMatchedClusters3D(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::Cluster *const pCluster3) const = 0;

    /**
     *  @brief  Calculate Pfo properties from proto particle
     *
     *  @param  protoParticle the input proto particle
     *  @param  pfoParameters the output pfo parameters
     */
    virtual void SetPfoParameters(const CosmicRayBaseMatchingAlgorithm::Particle &protoParticle,
        PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const = 0;

private:
    /**
     *  @brief Get a vector of available clusters
     *
     *  @param inputClusterListName the input name of the cluster list
     *  @param clusterVector the output vector of available clusters
     */
    pandora::StatusCode GetAvailableClusters(const std::string inputClusterListName, pandora::ClusterVector &clusterVector) const;
    /**
     *  @brief Match sets of clusters from two views
     *
     *  @param clusterVector1 the vector of clusters from the first view
     *  @param clusterVector2 the vector of clusters from the second view
     *  @param matchedClusters12 the map of cluster matches
     */
    void MatchClusters(const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2,
        CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters12) const;

    /**
     *  @brief Match clusters from three views and form into particles
     *
     *  @param matchedClusters12 the map of matches between the view 1 and view 2
     *  @param matchedClusters23 the map of matches between the view 2 and view 3
     *  @param matchedClusters31 the map of matches between the view 3 and view 1
     *  @param particleList the output list of particles
     */
    void MatchThreeViews(const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
        const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters23,
        const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters31,
        ParticleList &particleList) const;

    /**
     *  @brief Match clusters from two views and form into particles
     *
     *  @param matchedClusters12 the map of matches between the view 1 and view 2
     *  @param matchedClusters23 the map of matches between the view 2 and view 3
     *  @param matchedClusters31 the map of matches between the view 3 and view 1
     *  @param particleList the output list of particles
     */
    void MatchTwoViews(const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
        const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters23,
        const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters31,
        ParticleList &particleList) const;

    /**
     *  @brief Match clusters from two views and form into particles
     *
     *  @param matchedClusters12 the map of matches between the first and second views
     *  @param particleList the output list of particles
     */
    void MatchTwoViews(const CosmicRayBaseMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
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

    std::string    m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string    m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string    m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string    m_outputPfoListName;            ///< The name of the output PFO list
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_BASE_MATCHING_ALGORITHM_H
