/**
 *  @file   LArContent/include/LArClusterAssociation/ClusterExtensionAlgorithm.h
 * 
 *  @brief  Header file for the cluster extension algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_EXTENSION_ALGORITHM_H
#define LAR_CLUSTER_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  ClusterExtensionAlgorithm class
 */
class ClusterExtensionAlgorithm : public pandora::Algorithm
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
    void GetListOfCleanClusters( const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector );
    
    /**
     *  @brief  Form associations between pointing clusters
     * 
     *  @param  clusterI the first pointing cluster
     *  @param  clusterJ the second pointing cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix( const LArPointingCluster& clusterI, const LArPointingCluster& clusterJ, ClusterAssociationMatrix& clusterAssociationMatrix);

    /**
     *  @brief  Form associations between pointing cluster vertices
     * 
     *  @param  clusterVertexI the first pointing cluster vertex
     *  @param  clusterVertexJ the second pointing cluster vertex
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix( const LArPointingCluster::Vertex& clusterVertexI, const LArPointingCluster::Vertex& clusterVertexJ, ClusterAssociationMatrix& clusterAssociationMatrix );

    /**
     *  @brief  Form associations between pointing cluster vertices
     * 
     *  @param  clusterVertexI the first pointing cluster vertex
     *  @param  clusterVertexJ the second pointing cluster vertex
     *  @param  clusterAssociation the cluster association information
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix( const LArPointingCluster::Vertex& clusterVertexI, const LArPointingCluster::Vertex& clusterVertexJ, ClusterAssociation thisAssociation, ClusterAssociationMatrix& clusterAssociationMatrix );

    /**
     *  @brief  Determine whether two clusters are associated and should be merged
     * 
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     *
     *  @return boolean
     */
    bool AreClustersAssociated( pandora::Cluster* pCluster1, pandora::Cluster* pCluster2, ClusterAssociationMatrix& clusterAssociationMatrix );


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

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterExtensionAlgorithm::ClusterAssociation::ClusterAssociation() : 
    m_parent(NONE), m_daughter(NONE), m_association(NOTHING), m_strength(UNASSOCIATED), m_fom(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterExtensionAlgorithm::ClusterAssociation::ClusterAssociation( ClusterExtensionAlgorithm::ClusterAssociation::VertexType parent, ClusterExtensionAlgorithm::ClusterAssociation::VertexType daughter, ClusterExtensionAlgorithm::ClusterAssociation::AssociationType association, ClusterExtensionAlgorithm::ClusterAssociation::StrengthType strength, float fom) :
    m_parent(parent), m_daughter(daughter), m_association(association), m_strength(strength), m_fom(fom)
{
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

ClusterExtensionAlgorithm::ClusterAssociation::VertexType ClusterExtensionAlgorithm::ClusterAssociation::GetParent()
{
    return m_parent;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------
   
ClusterExtensionAlgorithm::ClusterAssociation::VertexType ClusterExtensionAlgorithm::ClusterAssociation::GetDaughter()
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterExtensionAlgorithm::ClusterAssociation::AssociationType ClusterExtensionAlgorithm::ClusterAssociation::GetAssociation()
{
    return m_association;
}
   
//------------------------------------------------------------------------------------------------------------------------------------------
      
ClusterExtensionAlgorithm::ClusterAssociation::StrengthType ClusterExtensionAlgorithm::ClusterAssociation::GetStrength()
{
    return m_strength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterExtensionAlgorithm::ClusterAssociation::GetFigureOfMerit()
{
    return m_fom;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_EXTENSION_ALGORITHM_H
