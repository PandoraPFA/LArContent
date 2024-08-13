/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/VertexBasedPfoRecoveryAlgorithm.h
 *
 *  @brief  Header file for the vertex-based particle recovery algorithm
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_BASED_PFO_RECOVERY_ALGORITHM_H
#define LAR_VERTEX_BASED_PFO_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  VertexBasedPfoRecoveryAlgorithm class
 */
class VertexBasedPfoRecoveryAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexBasedPfoRecoveryAlgorithm();

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
        Particle(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

        const pandora::Cluster *m_pClusterU; ///< Address of cluster in U view
        const pandora::Cluster *m_pClusterV; ///< Address of cluster in V view
        const pandora::Cluster *m_pClusterW; ///< Address of cluster in W view
    };

    typedef std::vector<Particle> ParticleList;

    /**
     *  @brief Get a vector of available clusters
     *
     *  @param inputClusterListName the input vector of the cluster list names
     *  @param clusterVector the output vector of available clusters
     */
    pandora::StatusCode GetAvailableClusters(const pandora::StringVector inputClusterListName, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Build the map of sliding fit results
     *
     *  @param  clusterVector the vector of selected clusters
     *  @param  halfWindowLayers the half-window to use for the sliding fits
     *  @param  slidingFitResultMap the sliding fit result map
     */
    void BuildSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief  Select clusters in proximity to reconstructed vertex
     *
     *  @param  pVertex  the input vertex
     *  @param  slidingFitResultMap  the mapping between clusters and sliding fit results
     *  @param  inputClusters  the input vector of clusters
     *  @param  outputClusters  the output vector of clusters
     */
    void SelectVertexClusters(const pandora::Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
        const pandora::ClusterVector &inputClusters, pandora::ClusterVector &outputClusters) const;

    /**
     *  @brief  Match clusters from three views
     *
     *  @param  pVertex  the input vertex
     *  @param  slidingFitResultMap  the mapping between clusters and sliding fit results
     *  @param  selectedClusters  the input vertex clusters
     *  @param  vetoList  the list of matched clusters
     *  @param  particleList the output list of matched clusters
     */
    void MatchThreeViews(const pandora::Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
        const pandora::ClusterVector &selectedClusters, pandora::ClusterSet &vetoList, ParticleList &particleList) const;

    /**
     *  @brief  Match clusters from two views
     *
     *  @param  pVertex  the input vertex
     *  @param  slidingFitResultMap  the mapping between clusters and sliding fit results
     *  @param  selectedClusters  the input vertex clusters
     *  @param  vetoList  the list of matched clusters
     *  @param  particleList  the output list of matched clusters
     */
    void MatchTwoViews(const pandora::Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
        const pandora::ClusterVector &selectedClusters, pandora::ClusterSet &vetoList, ParticleList &particleList) const;

    /**
     *  @brief  Get best-matched triplet of clusters from a set of input cluster vectors
     *
     *  @param  pVertex  the input vertex
     *  @param  slidingFitResultMap  the mapping between clusters and sliding fit results
     *  @param  clusters1  the clusters in the first view
     *  @param  clusters2  the clusters in the second view
     *  @param  clusters3  the clusters in the third view
     *  @param  pBestCluster1  the best-matched cluster from the first view
     *  @param  pBestCluster2  the best-matched cluster from the second view
     *  @param  pBestCluster3  the best-matched cluster from the third view
     *  @param  chi2  the chi-squared metric from the best match
     */
    void GetBestChi2(const pandora::Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
        const pandora::ClusterVector &clusters1, const pandora::ClusterVector &clusters2, const pandora::ClusterVector &clusters3,
        const pandora::Cluster *&pBestCluster1, const pandora::Cluster *&pBestCluster2, const pandora::Cluster *&pBestCluster3, float &chi2) const;

    /**
     *  @brief  Get best-matched pair of clusters from a set of input cluster vectors
     *
     *  @param  pVertex  the input vertex
     *  @param  slidingFitResultMap  the mapping between clusters and sliding fit results
     *  @param  clusters1 the clusters in the first view
     *  @param  clusters2 the clusters in the second view
     *  @param  pBestCluster1 the best-matched cluster from the first view
     *  @param  pBestCluster2 the best-matched cluster from the second view
     *  @param  chi2 the chi-squared metric from the best match
     */
    void GetBestChi2(const pandora::Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap, const pandora::ClusterVector &clusters1,
        const pandora::ClusterVector &clusters2, const pandora::Cluster *&pBestCluster1, const pandora::Cluster *&pBestCluster2, float &chi2) const;

    /**
     *  @brief  Merge two pointing clusters and return chi-squared metric giving consistency of matching
     *
     *  @param  pVertex  the input vertex
     *  @param  pointingCluster1  the first pointing cluster
     *  @param  pointingCluster2  the second pointing cluster
     */
    float GetChi2(const pandora::Vertex *const pVertex, const LArPointingCluster &pointingCluster1, const LArPointingCluster &pointingCluster2) const;

    /**
     *  @brief  Merge three clusters between views and return chi-squared metric giving consistency of matching
     *
     *  @param  pVertex  the input vertex
     *  @param  pointingCluster1  the first pointing cluster
     *  @param  pointingCluster2  the second pointing cluster
     *  @param  pointingCluster3  the third pointing cluster
     */
    float GetChi2(const pandora::Vertex *const pVertex, const LArPointingCluster &pointingCluster1,
        const LArPointingCluster &pointingCluster2, const LArPointingCluster &pointingCluster3) const;

    /**
     *  @brief  Select cluster which haven't been vetoed
     *
     *  @param  vetoList  the list of vetoed clusters
     *  @param  inputVector  the input vector of clusters
     *  @param  outputVector  the output vector of clusters
     */
    void SelectAvailableClusters(const pandora::ClusterSet &vetoList, const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;

    /**
     *  @brief  Select clusters of a specified hit type
     *
     *  @param  hitType  the specified hit type
     *  @param  inputVector  the input vector of clusters
     *  @param  outputVector  the output vector of clusters
     */
    void SelectClusters(const pandora::HitType hitType, const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;

    /**
     *  @brief  Find nearest end of pointing cluster to a specified position vector
     *
     *  @param  vertex the input position
     *  @param  cluster the input cluster
     */
    const LArPointingCluster::Vertex &GetInnerVertex(const pandora::CartesianVector &vertex, const LArPointingCluster &cluster) const;

    /**
     *  @brief  Find furthest end of pointing cluster from a specified position vector
     *
     *  @param  vertex the input position
     *  @param  cluster the input pointing cluster
     */
    const LArPointingCluster::Vertex &GetOuterVertex(const pandora::CartesianVector &vertex, const LArPointingCluster &cluster) const;

    /**
     *  @brief  Build particle flow objects from matched clusters
     *
     *  @param  particleList the input list of matched clusters
     */
    void BuildParticles(const ParticleList &particleList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputClusterListNames; ///< The list of input cluster list names
    std::string m_outputPfoListName;               ///< The name of the output pfo list

    unsigned int m_slidingFitHalfWindow; ///<
    float m_maxLongitudinalDisplacement; ///<
    float m_maxTransverseDisplacement;   ///<
    float m_twoViewChi2Cut;              ///<
    float m_threeViewChi2Cut;            ///<
    bool m_lowEnergyWorkflow;            ///<
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_BASED_PFO_RECOVERY_ALGORITHM_H
