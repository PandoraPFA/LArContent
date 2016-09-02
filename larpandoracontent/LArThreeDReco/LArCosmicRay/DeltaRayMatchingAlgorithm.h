/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h
 *
 *  @brief  Header file for the delta ray matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

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

    /**
     *  @brief  Default constructor
     */
    DeltaRayMatchingAlgorithm();

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
         *  @param  pCluster1 the first cluster
         *  @param  pCluster2 the second cluster
         *  @param  pCluster3 the third cluster
         *  @param  pPfo the parent Pfo
         */
        Particle(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
            const pandora::ParticleFlowObject *const pPfo);

        /**
         *  @brief  Get cluster in U view
         */
        const pandora::Cluster *GetClusterU() const;

        /**
         *  @brief  Get cluster in V view
         */
        const pandora::Cluster *GetClusterV() const;

        /**
         *  @brief  Get cluster in W view
         */
        const pandora::Cluster *GetClusterW() const;

        /**
         *  @brief  Get parent Pfo
         */
        const pandora::ParticleFlowObject *GetParentPfo() const;

        /**
         *  @brief  Get number of views
         */
        unsigned int GetNViews() const;

        /**
         *  @brief  Get number of calo hits
         */
        unsigned int GetNCaloHits() const;

        /**
         *  @brief  Get the total energy
         */
        float GetLengthSquared() const;

    private:
        const pandora::Cluster             *m_pClusterU;    ///< Address of cluster in U view
        const pandora::Cluster             *m_pClusterV;    ///< Address of cluster in V view
        const pandora::Cluster             *m_pClusterW;    ///< Address of cluster in W view
        const pandora::ParticleFlowObject  *m_pParentPfo;   ///< Address of parent Pfo
    };

    typedef std::vector<Particle> ParticleList;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::unordered_map<const pandora::Cluster*, pandora::ClusterList> ClusterToClustersMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;

    /**
     *  @brief  Initialize nearby cluster maps
     */
    void InitializeNearbyClusterMaps();

    /**
     *  @brief  Initialize a nearby cluster map with details relating to a specific cluster list
     * 
     *  @param  clusterListName the cluster list name
     *  @param  nearbyClustersMap to receive the nearby clusters map
     */
    void InitializeNearbyClusterMap(const std::string &clusterListName, ClusterToClustersMap &nearbyClusters);

    /**
     *  @brief  Clear nearby cluster maps
     */
    void ClearNearbyClusterMaps();

    /**
     *  @brief  Get a vector of all Pfos in the provided input Pfo lists
     *
     *  @param  inputPfoListName the input Pfo list name
     *  @param  pfoVector the output vector of Pfos
     */
    void GetAllPfos(const std::string &inputPfoListName, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Get a vector of track-like Pfos in the provided input Pfo lists
     *
     *  @param  inputPfoListName the input Pfo list name
     *  @param  pfoVector the output vector of Pfos
     */
    void GetTrackPfos(const std::string &inputPfoListName, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Get a vector containing all available input clusters in the provided cluster list, storing sliding linear fits
     *          in the algorithm cache.
     *
     *  @param  clusterListName the vector of cluster list names
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetClusters(const std::string &clusterListName, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Match clusters using all three views
     */
    void ThreeViewMatching() const;

    /**
     *  @brief  Match clusters using pairs of views
     */
    void TwoViewMatching() const;

    /**
     *  @brief  Match clusters using single views
     */
    void OneViewMatching() const;

    /**
     *  @brief  Match clusters using all three views
     *
     *  @param  clusters1 the list of clusters in the first view
     *  @param  clusters2 the list of clusters in the second view
     *  @param  clusters3 the list of clusters in the third view
     *  @param  particleList the output list of particles
     */
    void ThreeViewMatching(const pandora::ClusterVector &clusters1, const pandora::ClusterVector &clusters2,
        const pandora::ClusterVector &clusters3, ParticleList &particleList) const;

    /**
     *  @brief  Match clusters using a pair of views
     *
     *  @param  clusters1 the list of clusters in the first view
     *  @param  clusters2 the list of clusters in the second view
     *  @param  particleList the output list of particles
     */
    void TwoViewMatching(const pandora::ClusterVector &clusters1, const pandora::ClusterVector &clusters2,
        ParticleList &particleList) const;

    /**
     *  @brief  Match clusters using a single view
     *
     *  @param  clusters the list of clusters in the provided view
     *  @param  particleList the output list of particles
     */
    void OneViewMatching(const pandora::ClusterVector &clusters, ParticleList &particleList) const;

    /**
     *  @brief Resolve any ambiguities between candidate particles
     *
     *  @param inputParticles the input list of candidate particles
     *  @param outputParticles the output list of candidate particles
     */
    void SelectParticles(const ParticleList &inputParticles, ParticleList &outputParticles) const;

    /**
     *  @brief Build new particle flow objects
     *
     *  @param particleList the list of candidate particles
     */
    void CreateParticles(const ParticleList &particleList) const;

    /**
     *  @brief  Find best Pfo to associate a UVW triplet
     *
     *  @param  pointer to U view cluster
     *  @param  pointer to V view cluster
     *  @param  pointer to W view cluster
     *  @param  pointer to best Pfo
     */
    void FindBestParentPfo(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        const pandora::ParticleFlowObject *&pBestPfo) const;

    /**
     *  @brief  Look at consistency of a combination of clusters
     *
     *  @param  pCluster1 pointer to first luster
     *  @param  pCluster2 pointer to second cluster
     *  @param  pCluster3 pointer to third cluster
     */
    bool AreClustersMatched(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::Cluster *const pCluster3) const;

    /**
     *  @brief Get displacementr between cluster and particle flow object
     *
     *  @param pCluster pointer to cluster
     *  @param pPfo pointer to particle flow object
     */
    float GetDistanceSquaredToPfo(const pandora::Cluster *const pCluster, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Create a new Pfo from an input cluster list and set up a parent/daughter relationship
     *
     *  @param  clusterList the list of clusters
     *  @param  pParentPfo address of the parent pfo
     */
    void CreateDaughterPfo(const pandora::ClusterList &clusterList, const pandora::ParticleFlowObject *const pParentPfo) const;
  
    /**
     *  @brief  Merge an input cluster list with an existing daughter Pfo
     *
     *  @param  clusterList the list of clusters
     *  @param  pParentPfo address of the parent pfo
     */
    void AddToDaughterPfo(const pandora::ClusterList &clusterList, const pandora::ParticleFlowObject *const pParentPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_parentPfoListName;          ///< The parent pfo list name
    std::string             m_daughterPfoListName;        ///< The daughter pfo list name for new daughter particles

    std::string             m_inputClusterListNameU;      ///< The input cluster list name for the u view
    std::string             m_inputClusterListNameV;      ///< The input cluster list name for the v view
    std::string             m_inputClusterListNameW;      ///< The input cluster list name for the w view

    unsigned int            m_minCaloHitsPerCluster;      ///< The min number of calo hits per candidate cluster
    float                   m_xOverlapWindow;             ///< The maximum allowed displacement in x position
    float                   m_distanceForMatching;        ///< The maximum allowed distance between tracks and delta rays
    float                   m_pseudoChi2Cut;              ///< Pseudo chi2 cut for three view matching

    float                   m_searchRegion1D;             ///< Search region, applied to each dimension, for look-up from kd-trees
    ClusterToClustersMap    m_nearbyClustersU;            ///< The nearby clusters map for the u view
    ClusterToClustersMap    m_nearbyClustersV;            ///< The nearby clusters map for the v view
    ClusterToClustersMap    m_nearbyClustersW;            ///< The nearby clusters map for the w view
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *DeltaRayMatchingAlgorithm::Particle::GetClusterU() const
{
    return m_pClusterU;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline const pandora::Cluster *DeltaRayMatchingAlgorithm::Particle::GetClusterV() const
{
    return m_pClusterV;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline const pandora::Cluster *DeltaRayMatchingAlgorithm::Particle::GetClusterW() const
{
    return m_pClusterW;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline const pandora::ParticleFlowObject *DeltaRayMatchingAlgorithm::Particle::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayMatchingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_MATCHING_ALGORITHM_H
