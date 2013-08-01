/**
 *  @file   LArContent/include/LArClusterAssociation/TransverseAssociationAlgorithm.h
 * 
 *  @brief  Header file for the cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_ASSOCIATION_ALGORITHM_H
#define LAR_TRANSVERSE_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterAssociation/ClusterAssociationAlgorithm.h"

#include "Helpers/ClusterHelper.h"

namespace lar
{

/**
 *  @brief  TransverseAssociationAlgorithm class
 */
class TransverseAssociationAlgorithm : public ClusterAssociationAlgorithm
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
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  LArTransverseCluster class
     */
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
        const pandora::Cluster *GetCluster(unsigned int icluster = 0) const;

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
        pandora::CartesianVector    m_innerVertex;
        pandora::CartesianVector    m_outerVertex;
        pandora::CartesianVector    m_direction;
        pandora::ClusterVector      m_clusterVector;
        float                       m_rms;
    };

    typedef std::vector<LArTransverseCluster> LArTransverseClusterList;
    typedef std::map<const pandora::Cluster*, LArTransverseCluster> LArTransverseClusterMap;
    typedef std::map<pandora::Cluster*, pandora::ClusterList> LArClusterMergeMap; // TODO: Put this somewhere more central...

    /**
     *  @brief Separate clusters into transverse and longitudinal categories
     * 
     *  @param inputClusters input vector of clusters
     *  @param transverseClusters output vector of transverse or small clusters (small clusters are used to build new transverse clusters)
     *  @param longitudinalClusters output vector of longitudinal and large clusters
     */
    void GetTransverseClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &transverseClusters, pandora::ClusterVector &longitudinalClusters) const;

    /**
     *  @brief Select clean set of transverse clusters
     * 
     *  @param transverseClusters input vector of selected transverse clusters
     *  @param longitudinalClusters input vector of selected longitudinal clusters
     *  @param seedClusters output vector of transverse clusters selected as clean clusters
     *  @param nonSeedClusters output vector of transverse clusters not selected as clean clusters
     */
    void GetSeedClusters(const pandora::ClusterVector &transverseClusters, const pandora::ClusterVector &longitudinalClusters, pandora::ClusterVector &seedClusters, pandora::ClusterVector &nonSeedClusters) const;

    /**
     *  @brief Find the clusters that are transversely associated with a target cluster
     * 
     *  @param pCluster the target cluster
     *  @param inputClusters input vector of clean clusters for this event
     *  @param outputClusters output vector of clusters transversely associated with target cluster
     */
    void GetTransverseAssociatedClusters(const pandora::Cluster *const pCluster, const pandora::ClusterVector &inputClusters, pandora::ClusterVector &outputClusters) const;

    /**
     *  @brief Build a map of transverse cluster objects
     *
     *  @param seedClusters input vector of clean transverse clusters
     *  @param transverseClusters input vector of all transverse clusters
     *  @param transverseClusterMap output map of transverse cluster objects 
     */
    void FillTransverseClusterMap(const pandora::ClusterVector &seedClusters, const pandora::ClusterVector &transverseClusters, LArTransverseClusterMap &transverseClusterMap) const;

    /**
     *  @brief Build the map of forward/backward associations between clusters
     *
     *  @param transverseClusterMap input map of transverse cluster objects 
     *  @param clusterAssociationMap output map of forward/backward associations
     */
    void FillClusterAssociationMap(const LArTransverseClusterMap &transverseClusterMap, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Define whether two transverse cluster objects are associated
     *
     *  @param transCluster1 the inner cluster
     *  @param transCluster2 the outer cluster
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const LArTransverseCluster &transCluster1, const LArTransverseCluster &transCluster2) const;

    /**
     *  @brief Define whether a transverse cluster object is associated with the vertex position of a second transverse cluster object
     *
     *  @param theCluster the target cluster
     *  @param theVertex the vertex position
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const LArTransverseCluster &theCluster, const pandora::CartesianVector& theVertex) const;

    /**
     *  @brief Define whether clusters are transversely associated
     *
     *  @param pCluster1 the first cluster
     *  @param pCluster2 the second cluster
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief Get minimum and maximum X coordinates for a given cluster
     *
     *  @param pCluster the input cluster
     *  @param minX the minimum X position
     *  @param maxX the maximum X position
     */
    void GetExtremalCoordinatesX(const pandora::Cluster *const pCluster, float& minX, float& maxX) const;

    unsigned int   m_clusterCaloHits;               ///< 

    float          m_clusterWindow;                 ///< 
    float          m_clusterAngle;                  ///< 
    float          m_clusterCosAngle;               ///< 
    float          m_clusterTanAngle;               ///< 

    float          m_minCosRelativeAngle;           ///< 
    float          m_maxTransverseSeparation;       ///< 
    float          m_maxLongitudinalSeparation;     ///< 

    float          m_minTransverseLength;           ///< 
    unsigned int   m_minTransverseLayers;           ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const unsigned int TransverseAssociationAlgorithm::LArTransverseCluster::GetNumClusters() const
{
    return m_clusterVector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *TransverseAssociationAlgorithm::LArTransverseCluster::GetCluster(unsigned int icluster) const
{
    return m_clusterVector.at(icluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TransverseAssociationAlgorithm::LArTransverseCluster::GetInnerVertex() const
{
    return m_innerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TransverseAssociationAlgorithm::LArTransverseCluster::GetOuterVertex() const
{
    return m_outerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TransverseAssociationAlgorithm::LArTransverseCluster::GetDirection() const
{
    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const float TransverseAssociationAlgorithm::LArTransverseCluster::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TransverseAssociationAlgorithm::Factory::CreateAlgorithm() const
{
    return new TransverseAssociationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRANSVERSE_ASSOCIATION_ALGORITHM_H
