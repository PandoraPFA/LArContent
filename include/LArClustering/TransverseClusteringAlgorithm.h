/**
 *  @file   LArContent/include/Clustering/TransverseClusteringAlgorithm.h
 * 
 *  @brief  Header file for the cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_CLUSTERING_ALGORITHM_H
#define LAR_TRANSVERSE_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  TransverseClusteringAlgorithm class
 */
class TransverseClusteringAlgorithm : public pandora::Algorithm
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

 
    class LArTransverseCluster
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  
         */
        LArTransverseCluster(const pandora::Cluster *pSeedCluster, const pandora::ClusterVector &associatedClusters);

        /**
         *  @brief  Get the number of clusters
         * 
         *  @return the number of clusters
         */
        const unsigned int GetNumClusters() const;

        /**
         *  @brief  Get the address of the central cluster
         * 
         *  @return the address of the central cluster
         */
        const pandora::Cluster *GetCluster( unsigned int icluster = 0 ) const;

        /**
         *  @brief  Get the inner vertex position
         * 
         *  @return the inner vertex position
         */
        const pandora::CartesianVector &GetInnerVertex() const;

        /**
         *  @brief  Get the outer vertex position
         * 
         *  @return the outer vertex position
         */
        const pandora::CartesianVector &GetOuterVertex() const;

        /**
         *  @brief  Get the direction
         * 
         *  @return the direction
         */
        const pandora::CartesianVector &GetDirection() const;

        /**
         *  @brief  Get the RMS
         * 
         *  @return the RMS
         */
        const float GetRms() const;

    private:

        pandora::CartesianVector m_innerVertex;
        pandora::CartesianVector m_outerVertex;
        pandora::CartesianVector m_direction;

        pandora::ClusterVector   m_clusterVector;

        float m_rms;
    };

    typedef std::vector<LArTransverseCluster> LArTransverseClusterList;

    typedef std::map<const pandora::Cluster*, LArTransverseCluster> LArTransverseClusterMap;

    typedef std::map<pandora::Cluster*, bool> LArClusterVetoMap;

    // *** TODO, MOVE TO A MORE CENTRAL LOCATION ***
    typedef std::map<pandora::Cluster*, pandora::ClusterList> LArClusterMergeMap;
    

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    void GetSortedClusters( const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector );
    void GetTransverseClusters( const pandora::ClusterVector &inputClusters, pandora::ClusterVector &transverseClusters, pandora::ClusterVector &longitudinalClusters );
    void GetSeedClusters( const pandora::ClusterVector &transverseClusters, const pandora::ClusterVector &longitudinalClusters, pandora::ClusterVector &seedClusters, pandora::ClusterVector &nonSeedClusters );

    void FillTransverseClusterMap( const pandora::ClusterVector &seedClusters, const pandora::ClusterVector &transverseClusters, LArTransverseClusterMap &transverseClusterList ); 
    void GetTransverseAssociatedClusters( const pandora::Cluster* pCluster, const pandora::ClusterVector &inputClusters, pandora::ClusterVector &outputClusters ); 
    bool IsTransverseAssociated( const pandora::Cluster* pCluster1, const pandora::Cluster* pCluster2 );         


    void FillClusterMergeMap( const LArTransverseClusterMap &transverseClusterMap, LArClusterMergeMap &clusterMergeMap );
    bool IsTransverseAssociated( const LArTransverseCluster &transCluster1, const LArTransverseCluster &transCluster2 );
    bool IsTransverseAssociated( const LArTransverseCluster &theCluster, const pandora::CartesianVector& theVertex );

    void CollectAssociatedClusters( pandora::Cluster *pSeedCluster, pandora::Cluster *pCurrentCluster, LArClusterMergeMap &clusterMergeMap, LArClusterVetoMap &clusterVetoMap, pandora::ClusterList &associatedClusterList );


    
    

    
    



    unsigned int   m_clusterLayers;

    float          m_clusterWindow;
    float          m_clusterAngle;
    float          m_clusterCosAngle;
    float          m_clusterTanAngle; 

    float          m_minCosRelativeAngle;  
    float          m_maxTransverseSeparation;
    float          m_maxLongitudinalSeparation;

    float          m_minTransverseLength;
    unsigned int   m_minTransverseLayers;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const unsigned int TransverseClusteringAlgorithm::LArTransverseCluster::GetNumClusters() const
{
    return m_clusterVector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *TransverseClusteringAlgorithm::LArTransverseCluster::GetCluster( unsigned int icluster ) const
{
    return m_clusterVector.at(icluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TransverseClusteringAlgorithm::LArTransverseCluster::GetInnerVertex() const
{
    return m_innerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TransverseClusteringAlgorithm::LArTransverseCluster::GetOuterVertex() const
{
    return m_outerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TransverseClusteringAlgorithm::LArTransverseCluster::GetDirection() const
{
    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const float TransverseClusteringAlgorithm::LArTransverseCluster::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TransverseClusteringAlgorithm::Factory::CreateAlgorithm() const
{
    return new TransverseClusteringAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRANSVERSE_CLUSTERING_ALGORITHM_H
