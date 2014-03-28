/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/OvershootTracksTool.h
 * 
 *  @brief  Header file for the overshoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef OVERSHOOT_TRACKS_TOOL_H
#define OVERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  OvershootTracksTool class
 */
class OvershootTracksTool : public TensorManipulationTool
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
     *  @brief  Find overshoot tracks
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindOvershootTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Select elements representing possible components of a two particles, merged due to an overshoot in the clustering
     * 
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectOvershootElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
        const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Build proto particles, splitting clusters merged due to an overshoot in the clustering
     * 
     *  @param  iteratorList list of iterators to relevant elements
     *  @param  protoParticleVector to be populated with proto particles for subsequent pfo construction
     */
    void BuildProtoParticle(const IteratorList &iteratorList, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for overshoot identification
     * 
     *  @param  eIter the iterator to the tensor element
     *  @param  usedClusters the list of used clusters
     */
    bool PassesElementCuts(TensorType::ElementList::const_iterator eIter, const pandora::ClusterList &usedClusters) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key tensor element
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key tensor element
    float           m_minLongitudinalImpactParameter;   ///< The min longitudinal impact parameter for connecting accompanying clusters
    float           m_maxVertexXSeparation;             ///< The max separation between accompanying clusters vertex x positions to make split
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *OvershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new OvershootTracksTool();
}

} // namespace lar

#endif // #ifndef OVERSHOOT_TRACKS_TOOL_H
