/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/LongTracksTool.h
 * 
 *  @brief  Header file for the long tracks tool class.
 * 
 *  $Log: $
 */
#ifndef LONG_TRACKS_TOOL_H
#define LONG_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  LongTracksTool class
 */
class LongTracksTool : public TensorManipulationTool
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
     *  @brief  Whether a long element shares clusters with any other long elements
     * 
     *  @param  iIter specifies the long element under consideration
     *  @param  iteratorList list of iterators to other long elements
     * 
     *  @return boolean
     */
    static bool HasLongDirectConnections(IteratorList::const_iterator iIter, const IteratorList &iteratorList);

    /**
     *  @brief  Whether a long element is significantly longer that other elements with which it shares a cluster
     * 
     *  @param  iIter specifies the long element under consideration
     *  @param  elementList the full list of connected tensor elements
     *  @param  minMatchedSamplingPointRatio the min ratio between 1st and 2nd highest msps for simple ambiguity resolution
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     */
    static bool IsLongerThanDirectConnections(IteratorList::const_iterator iIter, const TensorType::ElementList &elementList,
        const unsigned int minMatchedSamplingPointRatio, const pandora::ClusterList &usedClusters);

    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find long tracks, hidden by simple ambiguities in the tensor
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindLongTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Select a list of long track-like elements from a set of connected tensor elements
     * 
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectLongElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters,
        IteratorList &iteratorList) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *LongTracksTool::Factory::CreateAlgorithmTool() const
{
    return new LongTracksTool();
}

} // namespace lar

#endif // #ifndef LONG_TRACKS_TOOL_H
