/**
 *  @file   LArContent/LArThreeDReco/LArCosmicRay/CosmicRayTrackRecoveryAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray longitudinal track recovery algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TRACK_RECOVERY_ALGORITHM_H
#define LAR_COSMIC_RAY_TRACK_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArContent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayTrackRecoveryAlgorithm class
 */
class CosmicRayTrackRecoveryAlgorithm : public pandora::Algorithm

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

    /**
     *  @brief  Default constructor
     */
    CosmicRayTrackRecoveryAlgorithm();

private:    
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        pandora::ClusterList m_clusterList;
    };

    typedef std::vector<Particle> ParticleList;
    typedef std::unordered_map<const pandora::Cluster*, pandora::ClusterList> ClusterAssociationMap;

    /**
     *  @brief Get a vector of available clusters
     *
     *  @param inputClusterListName the input name of the cluster list
     *  @param clusterVector the output vector of available clusters
     */
    pandora::StatusCode GetAvailableClusters(const std::string &inputClusterListName, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief Select a set of clusters judged to be clean
     *
     *  @param inputVector the input vector of all available clusters
     *  @param outputVector the output vector of clean clusters
     */
    void SelectCleanClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;

    /**
     *  @brief  Build the map of sliding fit results
     *
     *  @param  clusterVector the input cluster vector
     *  @param  slidingFitResultMap the output sliding fit result map
     */
    void BuildSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief Match a pair of cluster vectors and populate the cluster association map
     *
     *  @param clusterVector1 the input vector of clusters from the first view
     *  @param clusterVector2 the input vector of clusters from the second view
     *  @param slidingFitResultMap the input map of sliding linear fit results
     *  @param clusterAssociationMap the output map of cluster associations
     */
    void MatchViews(const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2, 
        const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Match a seed cluster with a list of target clusters and populate the cluster association map
     *
     *  @param pSeedCluster the input seed cluster
     *  @param targetClusters the input list of target clusters
     *  @param slidingFitResultMap the input map of sliding linear fit results
     *  @param clusterAssociationMap the output map of cluster associations
     */
    void MatchClusters(const pandora::Cluster* const pSeedCluster, const pandora::ClusterVector &targetClusters, 
        const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterAssociationMap &clusterAssociationMap) const;
   
    /**
     *  @brief  Create candidate particles using three primary clusters
     *
     *  @param clusterVectorU input vector of clusters from the U view
     *  @param clusterVectorV input vector of clusters from the V view
     *  @param clusterVectorW input vector of clusters from the W view
     *  @param clusterAssociationMapUV map of cluster associations between the U and V views
     *  @param clusterAssociationMapVW map of cluster associations between the V and W views
     *  @param clusterAssociationMapWU map of cluster associations between the W and U views
     *  @param particleList the output list of candidate particles
     */
    void MatchThreeViews(const pandora::ClusterVector &clusterVectorU, const pandora::ClusterVector &clusterVectorV,
        const pandora::ClusterVector &clusterVectorW, const ClusterAssociationMap &clusterAssociationMapUV, 
        const ClusterAssociationMap &clusterAssociationMapVW, const ClusterAssociationMap &clusterAssociationMapWU,
        ParticleList &particleList) const;

    /**
     *  @brief  Create candidate particles using two primary clusters and one pair of broken clusters
     *  
     *  @param clusterVectorU input vector of clusters from the U view
     *  @param clusterVectorV input vector of clusters from the V view
     *  @param clusterVectorW input vector of clusters from the W view
     *  @param clusterAssociationMapUV map of cluster associations between the U and V views
     *  @param clusterAssociationMapVW map of cluster associations between the V and W views
     *  @param clusterAssociationMapWU map of cluster associations between the W and U views
     *  @param particleList the output list of candidate particles
     */
    void MatchTwoViews(const pandora::ClusterVector &clusterVectorU,  const pandora::ClusterVector &clusterVectorV,
        const pandora::ClusterVector &clusterVectorW, const ClusterAssociationMap &clusterAssociationMapUV, 
        const ClusterAssociationMap &clusterAssociationMapVW, const ClusterAssociationMap &clusterAssociationMapWU,
        ParticleList &particleList) const;
 
    /**
     *  @brief  Create candidate particles using one primary cluster and one pair of broken clusters
     *
     *  @param clusterVectorU input vector of clusters from the U view
     *  @param clusterVectorV input vector of clusters from the V view
     *  @param clusterVectorW input vector of clusters from the W view
     *  @param clusterAssociationMapUV map of cluster associations between the U and V views
     *  @param clusterAssociationMapVW map of cluster associations between the V and W views
     *  @param clusterAssociationMapWU map of cluster associations between the W and U views
     *  @param particleList the output list of candidate particles
     */
    void MatchOneView(const pandora::ClusterVector &clusterVectorU,  const pandora::ClusterVector &clusterVectorV,
        const pandora::ClusterVector &clusterVectorW, const ClusterAssociationMap &clusterAssociationMapUV, 
        const ClusterAssociationMap &clusterAssociationMapVW, const ClusterAssociationMap &clusterAssociationMapWU,
        ParticleList &particleList) const;

    /**
     *  @brief Build the list of clusters already used to create particles
     *
     *  @param particleList the current list of candidate particles
     *  @param vetoList the list of clusters that belong to the candidate particles
     */
    void BuildVetoList(const ParticleList &particleList, pandora::ClusterList &vetoList) const;

    /**
     *  @brief Remove particles with duplicate clusters
     *
     *  @param inputParticleList the input list of candidate particles
     *  @param outputParticleList the input list of candidate particles with duplications removed
     */
    void RemoveAmbiguities(const ParticleList &inputParticleList, ParticleList &outputParticleList) const; 

    /**
     *  @brief Merge broken clusters into a single cluster
     *
     *  @param inputClusterList the input list of broken clusters 
     *  @param outputClusterList the output list of merged clusters
     */
    void MergeClusters(const pandora::ClusterList &inputClusterList, pandora::ClusterList &outputClusterList) const; 

    /**
     *  @brief Build particle flow objects
     *
     *  @param particleList the input list of candidate particles
     */
    void BuildParticles(const ParticleList &particleList);

    float          m_clusterMinLength;          ///<
    float          m_clusterMinSpanZ;           ///<
    float          m_clusterMinOverlapX;        ///<
    float          m_clusterMaxDeltaX;          ///<

    std::string    m_inputClusterListNameU;     ///<
    std::string    m_inputClusterListNameV;     ///<
    std::string    m_inputClusterListNameW;     ///<
    std::string    m_outputPfoListName;         ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackRecoveryAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackRecoveryAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_TRACK_RECOVERY_ALGORITHM_H
