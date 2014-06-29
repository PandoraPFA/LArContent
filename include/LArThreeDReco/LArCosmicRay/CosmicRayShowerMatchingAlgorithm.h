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

    

    typedef std::map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;

    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterAssociationMap;



    void MergeClusters();


    void MatchClusters();


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
     *  @brief  Generate a mapping between pfos and clusters
     *
     *  @param  pfoVector the input vector of Pfos
     *  @param  clusterToPfoMap the output mapping from clusters to pfos
     */
    void GetPfoClusterMap(const pandora::PfoVector &pfoVector, ClusterToPfoMap &clusterToPfoMap) const;
 


    void SelectPfos(const pandora::PfoVector &inputVector, pandora::PfoVector &outputVector) const;



 

    /**
     *  @brief  Fill list of candidate particles from input vector of Pfos
     *
     *  @param  pfoVector the input vector of Pfos
     *  @param  particleList the ouput list of candidate particles
     */
    void MatchViews(const pandora::PfoVector &pfoVector, ParticleList &particleList) const;

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
     *  @brief  Build map of cluster associations from list of candidate particles
     *
     *  @param  particleList the list of candidate particles
     *  @param  clusterAssociationMap the output map of cluster associations
     */
    void BuildAssociationMap(const ParticleList &particleList, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Build map of cluster associations from two clusters in the same view
     *
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  clusterAssociationMap the output map of cluster associations
     */
    void BuildAssociationMap(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, 
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Merge the associated clusters
     *
     *  @param  clusterToPfoMap the mapping between clusters and current Pfos
     *  @param  clusterAssociationMap the mapping between clusters to be merged
     */
    void MergeClusters(const ClusterToPfoMap &clusterToPfoMap, const ClusterAssociationMap &clusterAssociationMap) const;



    void CollectAssociatedClusters(const ClusterAssociationMap &clusterAssociationMap, 
        ClusterAssociationMap &clusterMergeMap) const;


    void CollectAssociatedClusters(pandora::Cluster *pSeedCluster, pandora::Cluster *pCurrentCluster, 
        const ClusterAssociationMap &clusterAssociationMap, const pandora::ClusterList &clusterVetoList, 
        pandora::ClusterList &associatedClusterList) const;
  

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    std::string  m_inputPfoListName;        ///< The input pfo list name
    std::string  m_inputClusterListNameU;     ///< The input cluster list name for the u view
    std::string  m_inputClusterListNameV;     ///< The input cluster list name for the v view
    std::string  m_inputClusterListNameW;     ///< The input cluster list name for the w view


    unsigned int m_minCaloHitsPerCluster;     ///< The min number of calo hits per candidate cluster

    float        m_xOverlapWindow;
    float        m_distanceForMatching;
    float        m_pseudoChi2Cut;

   
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
