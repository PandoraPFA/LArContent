/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h
 *
 *  @brief  Header file for the transverse association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_ASSOCIATION_ALGORITHM_H
#define LAR_TRANSVERSE_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Helpers/ClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

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

    /**
     *  @brief  Separate input clusters by length
     *
     *  @param  inputClusters the input vector of clusters
     *  @param  shortClusters the output vector of short clusters
     *  @param  transverseMediumClusters the output vector of transverse medium clusters
     *  @param  longitudinalMediumClusters the output vector of longitudinal medium clusters
     *  @param  longClusters the output vector of all long clusters
     */
    void SortInputClusters(const pandora::ClusterVector &inputClusters, pandora::ClusterVector &shortClusters,
        pandora::ClusterVector &transverseMediumClusters, pandora::ClusterVector &longitudinalMediumClusters,
        pandora::ClusterVector &longClusters) const;

    /**
     *  @brief  Form associations between two input lists of cluster
     *
     *  @param  firstVector the first input vector of clusters
     *  @param  secondVector the second input vector of clusters
     *  @param  firstAssociationMap the map of associations between first and second cluster vectors
     *  @param  secondAssociationMap the reversed map of associations between first and cluster vectors
     */
    void FillAssociationMap(const pandora::ClusterVector &firstVector, const pandora::ClusterVector &secondVector,
        ClusterAssociationMap &firstAssociationMap, ClusterAssociationMap &secondAssociationMap) const;

    /**
     *  @brief  Form a reduced set of associations between two input lists of clusters
     *
     *  @param  firstVector the first input vector of clusters
     *  @param  secondVector the second input vector of clusters
     *  @param  clusterAssociationMap the output map of associations between clusters
     */
    void FillReducedAssociationMap(const pandora::ClusterVector &firstVector, const pandora::ClusterVector &secondVector,
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Create transverse cluster objects, these are protoclusters with a direction and inner/outer vertices
     *
     *  @param  inputClusters the input vector of clusters
     *  @param  inputAssociationMap the map of associations between input clusters
     *  @param  transverseClusterList the output vector of transverse cluster objects
     */
    void FillTransverseClusterList(const  pandora::ClusterVector &inputClusters, const ClusterAssociationMap &inputAssociationMap,
        TransverseClusterList &transverseClusterList) const;

    /**
     *  @brief  Form associations between transverse cluster objects
     *
     *  @param  transverseClusterList the input vector of transverse cluster objects
     *  @param  transverseAssociationMap the external map of associations between clusters
     *  @param  clusterAssociationMap the output map of associations between clusters
     */
    void FillTransverseAssociationMap(const TransverseClusterList &transverseClusterList, const ClusterAssociationMap &transverseAssociationMap,
        ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Find the clusters that are transversely associated with a target cluster
     *
     *  @param  pCluster the target cluster
     *  @param  inputAssociationMap the map of associations between clusters
     *  @param  outputClusters the output vector of clusters transversely associated with target cluster
     */
    void GetAssociatedClusters(pandora::Cluster *const pCluster, const ClusterAssociationMap &inputAssociationMap,
        pandora::ClusterVector &associatedClusters) const;

    /**
     *  @brief  Determine whether clusters are association
     *
     *  @param  isForward whether the association is forwards or backwards
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *
     *  @return boolean
     */
    bool IsAssociated(const bool isForward, const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Determine whether two clusters are within the same cluster window
     *
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Determine whether two transverse clusters are associated
     *
     *  @param  pTransverseCluster1 the first transverse cluster
     *  @param  pTransverseCluster2 the second transverse cluster
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster1, const LArTransverseCluster *const pTransverseCluster2) const;

    /**
     *  @brief  Determine whether one transverse cluster is associated with the vertex from a second transverse cluster
     *
     *  @param  pTransverseCluster the target cluster
     *  @param  theVertex the vertex position
     *
     *  @return boolean
     */
    bool IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster, const pandora::CartesianVector &testPosition) const;

    /**
     *  @brief  Determine whether two clusters are overlapping
     *
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *
     *  @return boolean
     */
    bool IsOverlapping(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Calculate the overall span in X for a clusters
     *
     *  @param  pCluster the target cluster
     */
    float GetTransverseSpan(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Calculate the overall span in Z for a clusters
     *
     *  @param  pCluster the target cluster
     */
    float GetLongitudinalSpan(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Calculate the overall span in X for a set of clusters
     *
     *  @param  pCluster the target cluster
     *  @param  associatedClusters the vector of associated clusters
     */
    float GetTransverseSpan(const pandora::Cluster *const pCluster, const pandora::ClusterVector &associatedClusters) const;

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
     *  @brief  Get minimum and maximum X or Z coordinates for a given cluster
     *
     *  @param  pCluster the input cluster
     *  @param  useX calculate extermal coordinates for X (rather than Z)
     *  @param  minXZ the minimum X or Z position
     *  @param  maxXZ the maximum X or Z position
     */
    void GetExtremalCoordinatesXZ(const pandora::Cluster *const pCluster, const bool useX, float &minXZ, float &maxXZ) const;

    /**
     *  @brief  Get extremal 2D coordinates for a given cluster (ordered by X)
     *
     *  @param  pCluster the input cluster
     *  @param  innerCoordinate the inner coordinate
     *  @param  outerCoordinate the outer coordinate
     */
    void GetExtremalCoordinatesX(const pandora::Cluster *const pCluster, pandora::CartesianVector &innerCoordinate,
        pandora::CartesianVector &outerCoordinate) const;

    /**
     *  @brief Remove double-counting from association map
     *
     *  @param inputAssociationMap the inputted association map
     *  @param outputAssociationMap the outputted association map
     */
    void FillReducedAssociationMap(const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const;

    /**
     *  @brief Use one map to block associations from another map
     *
     *  @param firstAssociationMap the first association map
     *  @param secondAssociationMap the second association map
     *  @param secondAssociationMap the second association map reversed
     *  @param clusterAssociationMap the outputted association map
     */
    void FillReducedAssociationMap(const ClusterAssociationMap &firstAssociationMap, const ClusterAssociationMap &secondAssociationMap,
        const ClusterAssociationMap &secondAssociationMapSwapped, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Symmetrise an association map
     *
     *  @param inputAssociationMap the inputted association map
     *  @param outputAssociationMap the outputted association map
     */
    void FillSymmetricAssociationMap(const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const;

    /**
     *  @brief Symmetrise and then remove double-counting from an association map
     *
     *  @param inputAssociationMap the inputted association map
     *  @param outputAssociationMap the outputted association map
     */
    void FinalizeClusterAssociationMap(const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const;

    float          m_firstLengthCut;                   ///<
    float          m_secondLengthCut;                  ///<

    float          m_clusterWindow;                    ///<
    float          m_clusterAngle;                     ///<
    float          m_clusterCosAngle;                  ///<
    float          m_clusterTanAngle;                  ///<

    float          m_maxTransverseOverlap;             ///<
    float          m_maxLongitudinalOverlap;           ///<
    float          m_maxProjectedOverlap;              ///<

    float          m_transverseClusterMinCosTheta;     ///<
    float          m_transverseClusterMinLength;       ///<
    float          m_transverseClusterMaxDisplacement; ///<
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
