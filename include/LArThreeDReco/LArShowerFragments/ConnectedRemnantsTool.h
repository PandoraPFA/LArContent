/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerFragments/ConnectedRemnantsTool.h
 *
 *  @brief  Header file for the connected remnants tool class.
 *
 *  $Log: $
 */
#ifndef CONNECTED_REMNANTS_TOOL_H
#define CONNECTED_REMNANTS_TOOL_H 1

#include "LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ConnectedRemnantsTool class
 */
class ConnectedRemnantsTool : public RemnantTensorTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    ConnectedRemnantsTool();

    bool Run(ThreeDRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor);

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
     *  @param  connectedElements  the input list of connected elements
     *  @param  usedClusters  the list of clusters already analysed
     *  @param  clusterListU  the output list of clusters for the U view
     *  @param  clusterListV  the output list of clusters for the V view
     *  @param  clusterListW  the output list of clusters for the W view
     */
    void GetClusters(const TensorType::ElementList &connectedElements, const pandora::ClusterList &usedClusters,
        pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Fill map of clusters to be merged
     *
     *  @param  clusterList the input list of clusters
     *  @param  clusterMergeMap the output map of cluster merges
     */
    void FillMergeMap(const pandora::Cluster *const pCluster, const pandora::ClusterList &clusterList, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Check whether all clusters in the list are spatially connected
     *
     *  @param  clusterList  the input cluster list
     */
    bool IsConnected(const pandora::ClusterList &clusterList) const;

    float m_maxClusterSeparation;       ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ConnectedRemnantsTool::Factory::CreateAlgorithmTool() const
{
    return new ConnectedRemnantsTool();
}

} // namespace lar_content

#endif // #ifndef CONNECTED_REMNANTS_TOOL_H
