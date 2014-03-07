/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ClearTracksTool.h
 * 
 *  @brief  Header file for the clear tracks tool class.
 * 
 *  $Log: $
 */
#ifndef CLEAR_TRACKS_TOOL_H
#define CLEAR_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  ClearTracksTool class
 */
class ClearTracksTool : public TensorManipulationTool
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

private:
    pandora::StatusCode Run(const SlidingFitResultMap &slidingFitResultMap, TensorType &overlapTensor, ProtoParticleVector &protoParticleVector);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Ambiguity function for three-track case
     * 
     *  @param  clusterListU cluster list U
     *  @param  clusterListV cluster list V
     *  @param  clusterListW cluster list W
     *  @param  pClusterU to receive the address of the unambiguous U cluster
     *  @param  pClusterV to receive the address of the unambiguous V cluster
     *  @param  pClusterW to receive the address of the unambiguous W cluster
     * 
     *  @return boolean
     */
    static bool TrackTrackTrackAmbiguity(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV, const pandora::ClusterList &clusterListW,
        pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearTracksTool::Factory::CreateAlgorithmTool() const
{
    return new ClearTracksTool();
}

} // namespace lar

#endif // #ifndef CLEAR_TRACKS_TOOL_H
