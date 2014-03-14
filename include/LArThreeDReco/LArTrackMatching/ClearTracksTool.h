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
    pandora::StatusCode Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     * 
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     */
    void CreateThreeDParticles(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList) const;

    /**
     *  @brief  Classify elements of a cluster list as track-like or shower-like (uses IsMipTrack cluster flags, set by previous algs).
     * 
     *  @param  clusterList the input cluster list
     *  @param  trackClusterList to receive the list of track-like clusters
     *  @param  showerClusterList to receive the list of shower-like clusters
     */
    static void ClassifyClusters(const pandora::ClusterList &clusterList, pandora::ClusterList &trackClusterList, pandora::ClusterList &showerClusterList);

    static bool TrackTrackTrackAmbiguity(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW,  pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW);

    static bool TrackTrackShowerAmbiguity(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW,  pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW);

    static bool TrackShowerShowerAmbiguity(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW,  pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW);

    static bool ShowerShowerShowerAmbiguity(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW,  pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW);

    /**
     *  @brief  Resolve simple ambiguities in the tensor, cases where it is clear which element corresponds to primary track overlap
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  elementList to receive the list of tensor elements for particle creation
     */
    void ResolveSimpleAmbiguities(const TensorType &overlapTensor, TensorType::ElementList &elementList) const;

    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearTracksTool::Factory::CreateAlgorithmTool() const
{
    return new ClearTracksTool();
}

} // namespace lar

#endif // #ifndef CLEAR_TRACKS_TOOL_H
