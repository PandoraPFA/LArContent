/**
 *  @file   LArContent/include/LArClusterSeedAssociation/BoundedClusterMergingAlgorithm.h
 * 
 *  @brief  Header file for the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H
#define LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  BoundedCluster class
 */

class BoundedCluster 
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pCluster address of the cluster
     */
    BoundedCluster( const pandora::Cluster* pCluster = NULL);
    BoundedCluster( const BoundedCluster& rhs );
    ~BoundedCluster();

    /**
     *  @brief  Build a new bounding box from an inputted cluster
     * 
     *  @param  pCluster address of the cluster
     */
    void BuildBoundingBox( const pandora::Cluster* pCluster);
     
    /**
     *  @brief  Are the inner and outer layers of a cluster enclosed in the bounding box?
     * 
     *  @param  pCluster address of the cluster
     */
    unsigned int NumberOfEnclosedEnds( const pandora::Cluster* pCluster);

private:

    bool*  fBoxFlag;
    float* fBoxMinX;
    float* fBoxMaxX;
    float* fBoxMinY;
    float* fBoxMaxY;

    unsigned int fBoxLayers;
    unsigned int fNumLayers;
    unsigned int fMinLayer;
    unsigned int fMaxLayer;

};

/**
 *  @brief  BoundedClusterMergingAlgorithm class
 */

class BoundedClusterMergingAlgorithm : public pandora::Algorithm
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

    std::string         m_seedClusterListName;      ///< The seed cluster list name
    std::string         m_nonSeedClusterListName;   ///< The non seed cluster list name

    float               m_vertexVetoRadius;
    unsigned int        m_minClusterSize;

private:
    void ConstructBoundingBox( const pandora::Cluster* pCluster );
    unsigned int ApplyBoundingBox( const pandora::Cluster* pCluster );

    bool PassesVertexVeto( const pandora::Cluster* seedCluster, const pandora::Cluster* candCluster );
    bool PassesLengthVeto( const pandora::Cluster* seedCluster, const pandora::Cluster* candCluster );

    void PrepareMerges();
    void PrepareAssociations();

    void SetGoodAssociation( pandora::Cluster* pCandidateCluster );
    void SetBadAssociation( pandora::Cluster* pCandidateCluster );

    void MakeAssociation( pandora::Cluster* pSeedCluster, pandora::Cluster* pCandidateCluster, unsigned int numEnds );

    typedef std::map<pandora::Cluster*,pandora::Cluster*> ClusterAssociationMap;
    typedef std::map<pandora::Cluster*,pandora::ClusterList> ClusterMergeMap;

    ClusterAssociationMap fGoodAssociations;
    ClusterAssociationMap fBadAssociations;
    ClusterAssociationMap fDoubleAssociations;
    ClusterAssociationMap fSingleAssociations;
    ClusterAssociationMap fRepeatedDoubleAssociations;
    ClusterAssociationMap fRepeatedSingleAssociations;

    ClusterMergeMap fClusterMergeMap; 

    BoundedCluster fBoundingBox;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BoundedClusterMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BoundedClusterMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H
