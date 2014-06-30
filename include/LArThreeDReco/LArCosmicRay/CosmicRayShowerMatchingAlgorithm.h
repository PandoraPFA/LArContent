/**
 *  @file   LArContent/include/LArThreeDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h
 * 
 *  @brief  Header file for the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  CosmicRayShowerMatchingAlgorithm class
 */
class CosmicRayShowerMatchingAlgorithm : public pandora::Algorithm
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
     *  @brief Run cluster merging
     */
    void MergeClusters() const;

    /**
     *  @brief Run cluster matching (which builds new Pfos)
     */
    void MatchClusters() const;

    /**
     *  @brief  Get a vector of Pfos in the provided input Pfo lists
     *
     *  @param  inputPfoListName the input Pfo list name
     *  @param  pfoVector the output vector of Pfos
     */
    void GetPfos(const std::string inputPfoListName, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Get a vector of selected clusters in the provided input cluster lists
     *
     *  @param  clusterListName the input cluster list name
     *  @param  clusterVector to receive the output cluster vector
     */
    void GetClusters(const std::string inputClusterListName, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Select candidate showers and delta rays for matching (and don't select tracks)
     *
     *  @param  inputVector the input vector of all Pfos
     *  @param  outputVector the output vector of selected Pfos
     */
    void SelectPfos(const pandora::PfoVector &inputVector, pandora::PfoVector &outputVector) const;

    /**
     *  @brief Delete Pfos
     *
     *  @param  pfoVector the input vector of Pfos
     */
    void DeletePfos(const pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Merge clusters together based on results of matching
     *
     *  @param  particleList the input list of cluster matches
     */
    void MergeClusters(const ParticleList &particleList) const;

    /**
     *  @brief  Build Pfos based on matched clusters
     *
     *  @param  particleList the input list of cluster matches
     */
    void MatchClusters(const ParticleList &particleList) const;

    /**
     *  @brief  Fill list of candidate particles from input vectors of clusters
     *
     *  @param  clusterVectorU the input vector of clusters from the U view 
     *  @param  clusterVectorV the input vector of clusters from the V view 
     *  @param  clusterVectorW the input vector of clusters from the W view
     *  @param  particleList the ouput list of candidate particles
     */
    void MatchViews(const pandora::ClusterVector &clusterVectorU, const pandora::ClusterVector &clusterVectorV, 
        const pandora::ClusterVector &clusterVectorW, ParticleList &particleList) const;

   /**
     *  @brief  Match a UVW combination of clusters
     *
     *  @param  pClusterU the address of the cluster in the U view
     *  @param  pClusterV the address of the cluster in the V view
     *  @param  pClusterW the address of the cluster in the W view
     *
     *  @return boolean
     */
    bool MatchViews(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, 
        const pandora::Cluster *const pClusterW) const;

    /**
     *  @brief Identify good merges from list of matched clusters
     *
     *  @param particleList the input list of candidate particles
     *  @param clusterAssociationMap the output map of merges between clusters
     */   
    void IdentifyMerges(const ParticleList &particleList, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Build map of cluster associations from two clusters in the same view
     *
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  clusterAssociationMap the output map of cluster associations
     */
    void CheckAssociation(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, 
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Identify good matches from a list of possible matches
     *
     *  @param inputParticleList the input list of possible matches
     *  @param outputParticleList the output list of good matches 
     */   
    void IdentifyMatches(const ParticleList &inputParticleList, ParticleList &outputParticleList) const;

    /**
     *  @brief Perform cluster merging
     *
     *  @param clusterAssociationMap the input map of cluster merges
     */   
    void ApplyMerges(const ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Build new Pfos
     *
     *  @param particleList the input list of candidate particles 
     */   
    void ApplyMatches(const ParticleList &particleList) const;
 
    /**
     *  @brief Collect up cluster associations
     *
     *  @param clusterAssociationMap the input map of cluster associations
     *  @param clusterMergeMap the ouput map of cluster merges
     */   
    void CollectAssociatedClusters(const ClusterAssociationMap &clusterAssociationMap, 
        ClusterAssociationMap &clusterMergeMap) const;

    /**
     *  @brief Collect up cluster associations
     *
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  pCurrentCluster pointer to the current cluster
     *  @param  clusterAssociationMap the map of cluster associations
     *  @param  clusterVetoList the list of clusters that have already been merged
     *  @param  associatedClusterList the output list of associated clusters
     */   
    void CollectAssociatedClusters(pandora::Cluster *pSeedCluster, pandora::Cluster *pCurrentCluster, 
        const ClusterAssociationMap &clusterAssociationMap, const pandora::ClusterList &clusterVetoList, 
        pandora::ClusterList &associatedClusterList) const;
  
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_inputPfoListName;          ///< The input pfo list name
    std::string  m_inputClusterListNameU;     ///< The input cluster list name for the u view
    std::string  m_inputClusterListNameV;     ///< The input cluster list name for the v view
    std::string  m_inputClusterListNameW;     ///< The input cluster list name for the w view

    unsigned int m_minCaloHitsPerCluster;     ///< The min number of calo hits per candidate cluster

    float        m_xOverlapWindow;            ///< The overlap window in x coordinate when matching views
    float        m_distanceForMatching;       ///< Maximum distance between primary and secondary Pfos
    float        m_pseudoChi2Cut;             ///< Selection cut on pseudoChi2 when matching views
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
