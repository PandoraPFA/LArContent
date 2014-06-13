/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray splitting algorithm class.
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
     *  @brief Generate a map of sliding linear fit results from a vector of clusters
     *
     *  @param clusterVector the input vector of clusters
     *  @param slidingFitResultMap the output map of sliding linear fit results
     */
    void AddToSlidingFitCache(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitCache) const;




    void MatchTracks(const TwoDSlidingFitResultMap &slidingFitCache, 
        const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2, 
        CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12) const;



    bool MatchTracks(const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2) const;


   
    void MatchThreeViews(const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12, 
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters23,  
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters31,
        ParticleList &particleList) const;


    void MatchTwoViews(const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12, 
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters23,  
        const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters31,
        ParticleList &particleList) const;


    void MatchTwoViews(const CosmicRayTrackMatchingAlgorithm::ClusterAssociationMap &matchedClusters12,
        ParticleList &particleList) const;


    void ResolveAmbiguities(const ParticleList &inputList, ParticleList &outputList) const;
  

    void BuildParticles(const ParticleList &particleList);


    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string    m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string    m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string    m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string    m_outputPfoListName;            ///< The name of the output PFO list
  
    float          m_clusterMinLength;             ///< minimum length of clusters for this algorithm

    float          m_minXOverlap;                  ///< requirement on minimum X overlap for associated clusters
    float          m_minXOverlapFraction;          ///< requirement on minimum X overlap fraction for associated clusters
    unsigned int   m_halfWindowLayers;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
