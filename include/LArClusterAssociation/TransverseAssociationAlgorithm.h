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
         *  @param  pSeedCluster
         *  @param  associatedClusters
         */
        LArTransverseCluster(pandora::Cluster *pSeedCluster, const pandora::ClusterVector &associatedClusters);

        /**
         *  @brief  Constructor
         * 
         *  @return the address of the seed cluster
         */
        pandora::Cluster *GetSeedCluster() const;

        /**
         *  @brief  Get the associated cluster vector
         * 
         *  @return the associated cluster vector
         */
        const pandora::ClusterVector &GetAssociatedClusters() const;

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

    private:
        pandora::Cluster           *m_pSeedCluster;
        pandora::ClusterVector      m_associatedClusters;
        pandora::CartesianVector    m_innerVertex;
        pandora::CartesianVector    m_outerVertex;
        pandora::CartesianVector    m_direction;
    };

    typedef std::vector<LArTransverseCluster*> TransverseClusterList;
    typedef std::map<pandora::Cluster*, LArTransverseCluster*> TransverseClusterMap;
    typedef std::map<pandora::Cluster*, pandora::ClusterList> LArClusterMergeMap;

    /**
     *  @brief  Select building blocks for new transverse clusters
     * 
     *  @param  inputClusters input vector of clusters
     *  @param  transverseClusters output vector of building blocks for transverse clusters
     *  @param  longitudinalClusters output vector of established longitudinal clusters
     */
    void SeparateInputClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &transverseClusters, 
        pandora::ClusterVector &longitudinalClusters) const;

    /**
     *  @brief  Create transverse cluster objects, these are protoclusters with a direction and inner/outer vertices
     *
     *  @param  transverseClusters vector of building blocks for transverse clusters
     *  @param  longitudinalClusters vector of established longitudinal clusters
     *  @param  transverseClusterList output list of transverse cluster objects
     */
    void FillTransverseClusterList(const pandora::ClusterVector &transverseClusters, const pandora::ClusterVector &longitudinalClusters,
        TransverseClusterList &transverseClusterList) const;

    /**
     *  @brief  Use transverse cluster objects to define forward/backward associations
     *
     *  @param  transverseClusterList input list of transverse cluster objects 
     *  @param  forwardMergeMap
     *  @param  backwardMergeMap
     */
    void FillClusterMergeMaps(const TransverseClusterList &transverseClusterList, LArClusterMergeMap &forwardMergeMap,
        LArClusterMergeMap &backwardMergeMap) const;

    /**
     *  @brief  Use cluster merge maps to fill the final cluster association map
     *
     *  @param  forwardMergeMap
     *  @param  backwardMergeMap
     *  @param  clusterAssociationMap output map of forward/backward associations
     */
    void FillClusterAssociationMap(LArClusterMergeMap &forwardMergeMap, LArClusterMergeMap &backwardMergeMap,
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Find the clusters that are transversely associated with a target cluster
     * 
     *  @param  pCluster the target cluster
     *  @param  inputClusters input vector of clean clusters for this event
     *  @param  outputClusters output vector of clusters transversely associated with target cluster
     */
    void GetAssociatedClusters(pandora::Cluster *const pCluster, const pandora::ClusterVector &transverseClusters,
        const pandora::ClusterVector &longitudinalClusters, pandora::ClusterVector &associatedClusters) const;

    /**
     *  @brief  Define whether clusters are transversely associated
     *
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Define whether two transverse cluster objects are associated
     *
     *  @param  pTransverseCluster1 the inner cluster
     *  @param  pTransverseCluster2 the outer cluster
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster1, const LArTransverseCluster *const pTransverseCluster2) const;

    /**
     *  @brief  Define whether a transverse cluster object is associated with the vertex position of a second transverse cluster object
     *
     *  @param  pTransverseCluster the target cluster
     *  @param  theVertex the vertex position
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster, const pandora::CartesianVector &theVertex) const;

    /**
     *  @brief  Calculate the overall span in X for a set of clusters
     *
     *  @param  pCluster the target cluster
     *  @param  associatedClusters the vector of associated clusters
     *
     *  @return float
     */
    float GetTransverseLength(const pandora::Cluster *const pCluster, const pandora::ClusterVector &associatedClusters) const;

    /**
     *  @brief  Get minimum and maximum X or Z coordinates for a given cluster
     *
     *  @param  pCluster the input cluster
     *  @param  useX calculate extermal coordinates for X (rather than Z)
     *  @param  minXZ the minimum X or Z position
     *  @param  maxXZ the maximum X or Z position
     */
    void GetExtremalCoordinatesXZ(const pandora::Cluster *const pCluster, const bool useX, float &minXZ, float &maxXZ) const; 

    /**
     *  @brief  Get minimum and maximum X coordinates for a given cluster
     *
     *  @param  pCluster the input cluster
     *  @param  minX the minimum X position
     *  @param  maxX the maximum X position
     */
    void GetExtremalCoordinatesX(const pandora::Cluster *const pCluster, float &minX, float &maxX) const; 

    /**
     *  @brief  Get minimum and maximum Z coordinates for a given cluster
     *
     *  @param  pCluster the input cluster
     *  @param  minZ the minimum Z position
     *  @param  maxZ the maximum Z position
     */
    void GetExtremalCoordinatesZ(const pandora::Cluster *const pCluster, float &minZ, float &maxZ) const;  

    /**
     *  @brief  Get projected X coordinate of a cluster for a given Z coordinate
     *
     *  @param  pCluster the input cluster
     *  @param  inputZ the inputted Z position
     *  @param  outputZ the outputted Z projection
     */
    void GetProjectedCoordinateX(const pandora::Cluster *const pCluster, const float &inputZ, float &outputX) const;  

  

  


    


    float          m_clusterWindow;                 ///< 
    float          m_clusterAngle;                  ///< 
    float          m_clusterCosAngle;               ///< 
    float          m_clusterTanAngle;               ///< 

    float          m_minCosRelativeAngle;           ///< 
    float          m_maxTransverseSeparation;       ///< 

    float          m_minTransverseDisplacement;     ///< 
    float          m_maxLongitudinalDisplacement;   ///< 

    float          m_transverseClusterMaxCaloHits;  ///< 
    float          m_transverseClusterMaxLength;    ///< 
    float          m_longitudinalClusterMinLength;  ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *TransverseAssociationAlgorithm::LArTransverseCluster::GetSeedCluster() const
{
    return m_pSeedCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterVector &TransverseAssociationAlgorithm::LArTransverseCluster::GetAssociatedClusters() const
{
    return m_associatedClusters;
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
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TransverseAssociationAlgorithm::Factory::CreateAlgorithm() const
{
    return new TransverseAssociationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRANSVERSE_ASSOCIATION_ALGORITHM_H
