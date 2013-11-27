/**
 *  @file   LArContent/include/LArClusterAssociation/LongitudinalExtensionAlgorithm.h
 * 
 *  @brief  Header file for the cluster extension algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H
#define LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H 1

#include "LArClusterAssociation/ClusterMergingAlgorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  LongitudinalExtensionAlgorithm class
 */
class LongitudinalExtensionAlgorithm : public ClusterMergingAlgorithm
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
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const;

    /**
     *  @brief  LongitudinalAssociation class
     */
    class LongitudinalAssociation
    {
    public:
        /**
         *  @brief  Vertex enumeration
         */
        enum VertexType
        {
	    INNER = 1,
            OUTER = 2
        };

        /**
         *  @brief  Association enumeration
         */
        enum AssociationType
        {
            UNASSOCIATED = 0,
            BRANCH = 1,
            EMISSION = 2,
            NODE = 3
        };

        /**
         *  @brief  Strength enumeration
         */
        enum StrengthType
        {
            NONE = 0,
            WEAK = 1,
            STANDARD = 2,
            STRONG = 3
        };

        /**
         *  @brief  Constructor
         * 
         *  @param  parent
         *  @param  daughter
         *  @param  association
         *  @param  strength
         *  @param  fom
         */
        LongitudinalAssociation(const VertexType parent, const VertexType daughter, const AssociationType association,
            const StrengthType strength, const float fom);

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
         *  @brief  Get strength
         * 
         *  @return the strength
         */
        StrengthType GetStrength() const;

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
        StrengthType    m_strength;         ///< 
        float           m_fom;              ///< 
    };

    typedef std::map<const pandora::Cluster*, LongitudinalAssociation> LongitudinalAssociationMap;
    typedef std::map<const pandora::Cluster*, LongitudinalAssociationMap> LongitudinalAssociationMatrix;

    /**
     *  @brief  Fill the longitudinal association matrix
     * 
     *  @param  clusterVector the input vector of clusters
     *  @param  longitudinalAssociationMatrix the matrix of longitudinal associations
     */
    void FillAssociationMatrix(const pandora::ClusterVector &clusterVector, LongitudinalAssociationMatrix &longitudinalAssociationMatrix) const;

    /**
     *  @brief  Form associations between pointing clusters
     * 
     *  @param  clusterI the first pointing cluster
     *  @param  clusterJ the second pointing cluster
     *  @param  clusterAssociationMatrix the matrix of longitudinal associations
     */
    void FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, LongitudinalAssociationMatrix &clusterAssociationMatrix) const;

    /**
     *  @brief  Form associations between pointing cluster vertices
     * 
     *  @param  clusterVertexI the first pointing cluster vertex
     *  @param  clusterVertexJ the second pointing cluster vertex
     *  @param  clusterAssociationMatrix the matrix of longitudinal associations
     */
    void FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, const bool useInnerI, const bool useInnerJ, LongitudinalAssociationMatrix &clusterAssociationMatrix) const;

    
     /**
     *  @brief  Use association matrix to form associations between clusters
     * 
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *  @param  longitudinalAssociationMatrix the matrix of longitudinal associations
     */
    bool AreClustersAssociated(pandora::Cluster *pCluster1, pandora::Cluster *pCluster2, const LongitudinalAssociationMatrix &longitudinalAssociationMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float   m_spatialResolution;    ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LongitudinalExtensionAlgorithm::LongitudinalAssociation::LongitudinalAssociation(const VertexType parent, const VertexType daughter,
        const AssociationType association, const StrengthType strength, const float fom) :
    m_parent(parent),
    m_daughter(daughter),
    m_association(association),
    m_strength(strength),
    m_fom(fom)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LongitudinalExtensionAlgorithm::LongitudinalAssociation::VertexType LongitudinalExtensionAlgorithm::LongitudinalAssociation::GetParent() const
{
    return m_parent;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

inline LongitudinalExtensionAlgorithm::LongitudinalAssociation::VertexType LongitudinalExtensionAlgorithm::LongitudinalAssociation::GetDaughter() const
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LongitudinalExtensionAlgorithm::LongitudinalAssociation::AssociationType LongitudinalExtensionAlgorithm::LongitudinalAssociation::GetAssociation() const
{
    return m_association;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LongitudinalExtensionAlgorithm::LongitudinalAssociation::StrengthType LongitudinalExtensionAlgorithm::LongitudinalAssociation::GetStrength() const
{
    return m_strength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LongitudinalExtensionAlgorithm::LongitudinalAssociation::GetFigureOfMerit() const
{
    return m_fom;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LongitudinalExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new LongitudinalExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H
