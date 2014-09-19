/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h
 *
 *  @brief  Header file for the cluster extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_EXTENSION_ALGORITHM_H
#define LAR_CLUSTER_EXTENSION_ALGORITHM_H 1

#include "LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterExtensionAlgorithm class
 */
class ClusterExtensionAlgorithm : public ClusterMergingAlgorithm
{
protected:
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMatrix) const;

    /**
     *  @brief  ClusterAssociation class
     */
    class ClusterAssociation
    {
    public:
        /**
         *  @brief  Vertex enumeration
         */
        enum VertexType
        {
            UNDEFINED = 0,
            INNER     = 1,
            OUTER     = 2
        };

        /**
         *  @brief  Association enumeration
         */
        enum AssociationType
        {
            NONE   = 0,
            WEAK   = 1,
            STRONG = 2
        };

        /**
         *  @brief  Constructor
         *
         *  @param  parent
         *  @param  daughter
         *  @param  association
         *  @param  fom
         */
        ClusterAssociation(const VertexType parent, const VertexType daughter, const AssociationType association, const float fom);

        /**
         *  @brief  Get parent
         *
         *  @return the parent
         */
        VertexType GetParent() const;

        /**
         *  @brief  Get daughter
         *
         *  @return the daughter
         */
        VertexType GetDaughter() const;

        /**
         *  @brief  Get association
         *
         *  @return the association
         */
        AssociationType GetAssociation() const;

        /**
         *  @brief  Get figure of merit
         *
         *  @return the figure of merit
         */
        float GetFigureOfMerit() const;

    private:
        VertexType      m_parent;           ///<
        VertexType      m_daughter;         ///<
        AssociationType m_association;      ///<
        float           m_fom;              ///<
    };

    typedef std::map<const pandora::Cluster*, ClusterAssociation> ClusterAssociationMap;
    typedef std::map<const pandora::Cluster*, ClusterAssociationMap> ClusterAssociationMatrix;

    /**
     *  @brief  Fill the cluster association matrix
     *
     *  @param  clusterVector the input vector of clusters
     *  @param  clusterAssociationMatrix the matrix of associations
     */
    virtual void FillClusterAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const = 0;

     /**
     *  @brief  Fill the cluster merge map
     *
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     *  @param  clusterMergeMap the map of cluster merges
     */
    virtual void FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterExtensionAlgorithm::ClusterAssociation::ClusterAssociation(const VertexType parent, const VertexType daughter,
        const AssociationType association, const float fom) :
    m_parent(parent),
    m_daughter(daughter),
    m_association(association),
    m_fom(fom)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterExtensionAlgorithm::ClusterAssociation::VertexType ClusterExtensionAlgorithm::ClusterAssociation::GetParent() const
{
    return m_parent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterExtensionAlgorithm::ClusterAssociation::VertexType ClusterExtensionAlgorithm::ClusterAssociation::GetDaughter() const
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterExtensionAlgorithm::ClusterAssociation::AssociationType ClusterExtensionAlgorithm::ClusterAssociation::GetAssociation() const
{
    return m_association;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ClusterExtensionAlgorithm::ClusterAssociation::GetFigureOfMerit() const
{
    return m_fom;
}

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_EXTENSION_ALGORITHM_H
