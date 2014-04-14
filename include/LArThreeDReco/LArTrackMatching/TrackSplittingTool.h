/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/TrackSplittingTool.h
 * 
 *  @brief  Header file for the track splitting tool class.
 * 
 *  $Log: $
 */
#ifndef TRACK_SPLITTING_TOOL_H
#define TRACK_SPLITTING_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  TrackSplittingTool class
 */
class TrackSplittingTool : public TensorManipulationTool
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

    typedef std::map<pandora::Cluster*, pandora::CartesianPointList> SplitPositionMap;

    /**
     *  @brief  Find remaining tracks, hidden by spurious track segments (and maybe other ambiguities) in the tensor
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  splitPositionMap to receive the split position map
     */
    void FindTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, SplitPositionMap &splitPositionMap) const;

    typedef ThreeDTransverseTracksAlgorithm::IteratorList IteratorList;

    /**
     *  @brief  Select a list of the relevant elements from a set of connected tensor elements
     * 
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Whether a provided tensor element can be used to construct a pfo
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the tensor element
     *  @param  usedClusters the list of used clusters
     *  @param  splitPositionMap to receive the split position map
     */
    bool PassesChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element, pandora::ClusterList &usedClusters,
        SplitPositionMap &splitPositionMap) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (between long clusters and short cluster vs. shared overlap)
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TrackSplittingTool::Factory::CreateAlgorithmTool() const
{
    return new TrackSplittingTool();
}

} // namespace lar

#endif // #ifndef TRACK_SPLITTING_TOOL_H
