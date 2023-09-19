/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseSegmentMergingAlgorithm.h
 *
 *  @brief  Header file for the transverse segment merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_SEGMENT_MERGING_ALGORITHM_H
#define LAR_TRANSVERSE_SEGMENT_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseSegmentMergingAlgorithm class
 */
class TransverseSegmentMergingAlgorithm : public ClusterAssociationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TransverseSegmentMergingAlgorithm();

private:
    typedef std::map<const pandora::Cluster *, pandora::ClusterSet> ClusterToClustersMap;
    typedef std::map<const pandora::Cluster *, bool> ClusterUseMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Use extremal coordinates to identify nearby clusters for each cluster
     *
     *  @param  allClusters the list of all clusters
     *  @param  nearbyClusters the map of related clusters based on proximity
     */
    void GetNearbyClusterMap(const pandora::ClusterVector &allClusters, ClusterToClustersMap &nearbyClusters) const;

    /**
     *  @brief  Walks through the map of nearby clusters, starting from a given seed, to build chains of potential associations
     *
     *  @param  pSeedCluster the cluster with which to begin the chain of potential associations
     *  @param  nearbyClusters the input map of related clusters
     *  @param  usedClusters the map of clusters already used
     *  @param  associatedClusters the output set of clusters associated with the seed cluster
     */
    void GetIndirectAssociations(const pandora::Cluster *const pSeedCluster, const ClusterToClustersMap &nearbyClusters, ClusterUseMap &usedClusters,
        pandora::ClusterSet &associatedClusters) const;
};

} // namespace lar_content

#endif // #ifndef LAR_TRANSVERSE_SEGMENT_MERGING_ALGORITHM_H
