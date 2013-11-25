/**
 *  @file   LArContent/include/LArClusterSeedAssociation/ParallelClusterMergingAlgorithm.h
 * 
 *  @brief  Header file for the parallel cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PARALLEL_CLUSTER_MERGING_ALGORITHM_H
#define LAR_PARALLEL_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ParallelClusterMergingAlgorithm class
 */
class ParallelClusterMergingAlgorithm : public pandora::Algorithm
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
     *  @brief  ClusterProperties class
     */
    class ClusterProperties
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster address of the cluster
         */
        ClusterProperties(pandora::Cluster *const pCluster);

        /**
         *  @param  Get the cluster
         * 
         *  @return address of the cluster
         */
        pandora::Cluster *GetCluster() const;

        /**
         *  @param  Whether the cluster travels forwards in z
         * 
         *  @return boolean
         */
        bool IsForwardInZ() const;

        /**
         *  @param  Whether the cluster travels backwards in z
         * 
         *  @return boolean
         */
        bool IsBackwardInZ() const;

        /**
         *  @param  Whether the cluster is declared to be a muon
         * 
         *  @return boolean
         */
        bool IsMuon() const;

    private:
        pandora::Cluster   *m_pCluster;         ///< 
        bool                m_isForwardInZ;     ///< 
        bool                m_isBackwardInZ;    ///< 
        bool                m_isMuon;           ///< 
    };

    typedef std::vector<ClusterProperties> ClusterPropertiesList;

    /**
     *  @brief  Whether cluster can be declared as clean
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    bool IsCleanCluster(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Whether seed and target cluster overlap
     * 
     *  @param  pSeedCluster address of the seed cluster
     *  @param  pTargetCluster address of the target cluster
     * 
     *  @return boolean
     */
    bool IsParallel(const pandora::Cluster *const pSeedCluster, const pandora::Cluster *const pTargetCluster) const;

    /**
     *  @brief  Association checks on seed and target clusters
     * 
     *  @param  seedProperties the seed cluster properties
     *  @param  pTargetCluster address of the target cluster
     * 
     *  @return boolean
     */
    bool IsGoodMerge(const ClusterProperties &seedProperties, const pandora::Cluster *const pTargetCluster) const;

     /**
     *  @brief  Get closest distance squared between target cluster and seed vertex
     * 
     *  @param  seedProperties the seed cluster properties
     *  @param  pTargetCluster address of the target cluster
     */
    float GetClosestDistanceToVertexSquared(const ClusterProperties &seedProperties, const pandora::Cluster *const pTargetCluster) const;

    /**
     *  @brief  Get closest distance squared between target cluster and seed cluster
     * 
     *  @param  seedProperties the seed cluster properties
     *  @param  pTargetCluster address of the target cluster
     */
    float GetClosestDistanceSquared(const ClusterProperties &seedProperties, const pandora::Cluster *const pTargetCluster) const;

    /**
     *  @brief  Get closest distance squared between a cluster and a specified position
     * 
     *  @param  pCluster address of the cluster
     *  @param  position the position vector
     */
    float GetClosestDistanceSquared(const pandora::Cluster* pCluster, const pandora::CartesianVector &position) const;

    typedef std::map<pandora::Cluster*,pandora::ClusterList> ClusterMergeMap;

    std::string         m_seedClusterListName;                  ///< The seed cluster list name
    std::string         m_nonSeedClusterListName;               ///< The non seed cluster list name

    unsigned int        m_maxNumIterations;                     ///< 
    unsigned int        m_minClusterSize;                       ///< 
    float               m_clusterWindowLayers;                  ///< 

    float               m_clusterWindowRadius;                  ///< 
    float               m_clusterWindowCloseRadius;             ///< 
    float               m_vertexVetoRadius;                     ///< 

    float               m_clusterWindowRadiusSquared;           ///< 
    float               m_clusterWindowCloseRadiusSquared;      ///< 
    float               m_vertexVetoRadiusSquared;              ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParallelClusterMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParallelClusterMergingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ParallelClusterMergingAlgorithm::ClusterProperties::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ParallelClusterMergingAlgorithm::ClusterProperties::IsForwardInZ() const
{
    return m_isForwardInZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ParallelClusterMergingAlgorithm::ClusterProperties::IsBackwardInZ() const
{
    return m_isBackwardInZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ParallelClusterMergingAlgorithm::ClusterProperties::IsMuon() const
{
    return m_isMuon;
}

} // namespace lar

#endif // #ifndef LAR_PARALLEL_CLUSTER_MERGING_ALGORITHM_H
