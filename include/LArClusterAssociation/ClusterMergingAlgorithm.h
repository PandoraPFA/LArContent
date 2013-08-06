/**
 *  @file   LArContent/include/LArClusterAssociation/ClusterMergingAlgorithm.h
 * 
 *  @brief  Header file for the cluster extension algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
#define LAR_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ClusterMergingAlgorithm class
 */
class ClusterMergingAlgorithm : public pandora::Algorithm
{
protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

 
    /**
     *  ClusterAssociation class
     */

    class ClusterAssociation
    {
    public:   

        /**
         *  Vertex enumeration
         */
        enum VertexType
        {
            NONE = 0,
            INNER = 1,
            OUTER = 2
        };

        /**
         *  Association enumeration
         */
        enum AssociationType
        {
            NOTHING = 0,
            BRANCH = 1,
            EMISSION = 2,
            NODE = 3
        };

        /**
         *  Strength enumeration
         */
        enum StrengthType
        {
            UNASSOCIATED = 0,
            WEAK = 1,
            STANDARD = 2,
            STRONG = 3
        };

        ClusterAssociation();
        ClusterAssociation( StrengthType strength );
        ClusterAssociation( VertexType parent, VertexType daughter, 
                            AssociationType association, StrengthType strength, float fom );
 
        VertexType       GetParent();
        VertexType       GetDaughter();
        AssociationType  GetAssociation();
        StrengthType     GetStrength();
        float            GetFigureOfMerit();

    private:
        VertexType      m_parent;
        VertexType      m_daughter;
        AssociationType m_association;
        StrengthType    m_strength;   
        float           m_fom;
    };


    typedef std::map<const pandora::Cluster*, ClusterAssociation> ClusterAssociationMap;
    typedef std::map<const pandora::Cluster*, ClusterAssociationMap> ClusterAssociationMatrix;
  
    typedef std::map<pandora::Cluster*, bool> ClusterVetoMap;
    typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     * 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    virtual void GetListOfCleanClusters( const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector ) const = 0;
   

    /**
     *  @brief  Form associations between pointing clusters
     * 
     *  @param  clusterVector the vector of clean clusters
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    virtual void FillAssociationMatrix( const pandora::ClusterVector& clusterVector, ClusterAssociationMatrix& clusterAssociationMatrix) const = 0;


    /**
     *  @brief  Determine whether two clusters are associated and should be merged
     * 
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     *
     *  @return boolean
     */
    virtual bool AreClustersAssociated( pandora::Cluster* pCluster1, pandora::Cluster* pCluster2, ClusterAssociationMatrix& clusterAssociationMatrix ) const;

private:

  
   /**
     *  @brief  Collect up all clusters associations related to a given seed cluster
     * 
     *  @param  pSeedCluster pointer to the initial cluster
     *  @param  pCurrentCluster pointer to the current cluster 
     *  @param  clusterMergeMap the map of clusters to be merged
     *  @param  clusterVetoMap the map of clusters that have already been merged
     *  @param  associatedClusterList the list of associated clusters
     */
    void CollectAssociatedClusters( pandora::Cluster* pSeedCluster, pandora::Cluster* pCurrentCluster, ClusterMergeMap& clusterMergeMap, ClusterVetoMap& clusterVetoMap, pandora::ClusterList& associatedClusterList );
  
    
};



} // namespace lar

#endif // #ifndef LAR_CLUSTER_MERGING_ALGORITHM_H
