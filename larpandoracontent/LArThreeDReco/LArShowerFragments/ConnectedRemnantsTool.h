/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h
 *
 *  @brief  Header file for the connected remnants tool class.
 *
 *  $Log: $
 */
#ifndef CONNECTED_REMNANTS_TOOL_H
#define CONNECTED_REMNANTS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeViewRemnantsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ConnectedRemnantsTool class
 */
class ConnectedRemnantsTool : public RemnantTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ConnectedRemnantsTool();

    bool Run(ThreeViewRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identify candidate particles
     *
     *  @param  overlapTensor  the input overlap tensor
     *  @param  protoParticleVector  the output vector of candidate particles
     *  @param  clusterMergeMap  the output map of clusters to be merged
     */
    void FindConnectedShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Separate connected clusters into cluster lists by view
     *
     *  @param  connectedElements the input list of connected elements
     *  @param  usedClusters the list of clusters already analysed
     *  @param  clusterVectorU the output vector of clusters for the U view
     *  @param  clusterVectorV the output vector of clusters for the V view
     *  @param  clusterVectorW the output vector of clusters for the W view
     */
    void GetClusters(const TensorType::ElementList &connectedElements, const pandora::ClusterSet &usedClusters,
        pandora::ClusterVector &clusterVectorU, pandora::ClusterVector &clusterVectorV, pandora::ClusterVector &clusterVectorW) const;

    /**
     *  @brief  Fill map of clusters to be merged
     *
     *  @param  clusterVector the input vector of clusters
     *  @param  clusterMergeMap the output map of cluster merges
     */
    void FillMergeMap(const pandora::Cluster *const pCluster, const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Check whether all clusters in a list are spatially connected
     *
     *  @param  clusterVector the input cluster vector
     */
    bool IsConnected(const pandora::ClusterVector &clusterVector) const;

    float m_maxClusterSeparation; ///<
};

} // namespace lar_content

#endif // #ifndef CONNECTED_REMNANTS_TOOL_H
