/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/UndershootTracksTool.h
 * 
 *  @brief  Header file for the undershoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef UNDERSHOOT_TRACKS_TOOL_H
#define UNDERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  UndershootTracksTool class
 */
class UndershootTracksTool : public TensorManipulationTool
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
    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find undershoot tracks
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindUndershootTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Apply the changes cached in a proto particle vector and update the tensor accordingly
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  protoParticleVector the proto particle vector
     * 
     *  @return whether changes to the tensor have been made
     */
    bool ApplyChanges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Merge together clusters in a provided cluster list
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterList the cluster list
     * 
     *  @return whether changes to the tensor have been made
     */
    bool MakeClusterMerges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const pandora::ClusterList &clusterList) const;

    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Select elements representing possible components of a single particle, separated due to an undershoot in the clustering
     * 
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectUndershootElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
        const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Build a proto particle, merging together elements separated due to an undershoot in the clustering
     * 
     *  @param  iteratorList list of iterators to relevant elements
     *  @param  protoParticle to be populated with clusters for subsequent pfo construction
     */
    void BuildProtoParticle(const IteratorList &iteratorList, ProtoParticle &protoParticle) const;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for undershoot identification
     * 
     *  @param  eIter the iterator to the tensor element
     *  @param  usedClusters the list of used clusters
     */
    bool PassesElementCuts(TensorType::ElementList::const_iterator eIter, const pandora::ClusterList &usedClusters) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key tensor element
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key tensor element
    float           m_minLongitudinalImpactParameter;   ///< The minimum longitudinal impact parameter for connecting broken clusters
    float           m_maxTransverseImpactParameter;     ///< The maximum transverse impact parameter for connecting broken clusters
    float           m_minImpactParameterCosTheta;       ///< The minimum cos theta (angle between vertex directions) for connecting broken clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *UndershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new UndershootTracksTool();
}

} // namespace lar

#endif // #ifndef UNDERSHOOT_TRACKS_TOOL_H
