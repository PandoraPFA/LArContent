/**
 *  @file   LArContent/include/LArClusterAssociation/TrackExtensionAlgorithm.h
 * 
 *  @brief  Header file for the cluster extension algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_EXTENSION_ALGORITHM_H
#define LAR_TRACK_EXTENSION_ALGORITHM_H 1

#include "LArClusterAssociation/ClusterMergingAlgorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  TrackExtensionAlgorithm class
 */
class TrackExtensionAlgorithm : public ClusterMergingAlgorithm
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
    void GetListOfCleanClusters( const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector ) const;
    void FillAssociationMatrix( const pandora::ClusterVector& clusterVector, ClusterAssociationMatrix& clusterAssociationMatrix) const;
    bool AreClustersAssociated( pandora::Cluster* pCluster1, pandora::Cluster* pCluster2, ClusterAssociationMatrix& clusterAssociationMatrix ) const;


    /**
     *  @brief  Form associations between pointing clusters
     * 
     *  @param  clusterI the first pointing cluster
     *  @param  clusterJ the second pointing cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix( const LArPointingCluster& clusterI, const LArPointingCluster& clusterJ, ClusterAssociationMatrix& clusterAssociationMatrix) const;

    /**
     *  @brief  Form associations between pointing cluster vertices
     * 
     *  @param  clusterVertexI the first pointing cluster vertex
     *  @param  clusterVertexJ the second pointing cluster vertex
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix( const LArPointingCluster::Vertex& clusterVertexI, const LArPointingCluster::Vertex& clusterVertexJ, ClusterAssociationMatrix& clusterAssociationMatrix ) const;

  
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TrackExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new TrackExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRACK_EXTENSION_ALGORITHM_H
