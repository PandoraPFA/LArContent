/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h
 * 
 *  @brief  Header file for the track splitting tool class.
 * 
 *  $Log: $
 */
#ifndef TRACK_SPLITTING_TOOL_H
#define TRACK_SPLITTING_TOOL_H 1

#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TrackSplittingTool class
 */
class TrackSplittingTool : public TransverseTensorTool
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
    TrackSplittingTool();

    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  element the tensor element
         */
        Particle(const TensorType::Element &element);

        pandora::Cluster   *m_pLongCluster;             ///< Address of the long cluster
        pandora::Cluster   *m_pCluster1;                ///< Address of short cluster in view 1
        pandora::Cluster   *m_pCluster2;                ///< Address of short cluster in view 2
        float               m_longMinX;                 ///< The min x coordinate of the long cluster
        float               m_longMaxX;                 ///< The max x coordinate of the long cluster
        float               m_short1MinX;               ///< The min x coordinate of short cluster 1
        float               m_short1MaxX;               ///< The max x coordinate of short cluster 1
        float               m_short2MinX;               ///< The min x coordinate of short cluster 2
        float               m_short2MaxX;               ///< The max x coordinate of short cluster 2
    };

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find remaining tracks, hidden by spurious track segments (and maybe other ambiguities) in the tensor
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  splitPositionMap to receive the split position map
     */
    void FindTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, SplitPositionMap &splitPositionMap) const;

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
     *  @param  isMinX whether to look for track splits at min or max x coordinate
     *  @param  usedClusters the list of used clusters
     *  @param  splitPositionMap to receive the split position map
     */
    bool PassesChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element, const bool isMinX, pandora::ClusterList &usedClusters,
        SplitPositionMap &splitPositionMap) const;

    /**
     *  @brief  Check a candidate split position for consistency with the associated track cluster sliding linear fit
     * 
     *  @param  splitPosition the candidate split position
     *  @param  splitX the split x coordinate
     *  @param  longFitResult the sliding linear fit for the long cluster
     */
    bool CheckSplitPosition(const pandora::CartesianVector &splitPosition, const float splitX, const TwoDSlidingFitResult &longFitResult) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (between long clusters and short cluster vs. shared overlap)
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution

    float           m_maxShortDeltaXFraction;           ///< Max x distance between ends of two short clusters (measured as fraction of long cluster x length)
    float           m_maxAbsoluteShortDeltaX;           ///< Max x distance between ends of two short clusters (measured as an absolute distance)
    float           m_minLongDeltaXFraction;            ///< Min x distance between ends of short and long clusters (measured as fraction of long cluster x length)
    float           m_minAbsoluteLongDeltaX;            ///< Min x distance between ends of short and long clusters (measured as an absolute distance)
    float           m_minSplitToVertexProjection;       ///< Min projected distance between split position and either inner or outer vertex of long cluster
    float           m_maxSplitVsFitPositionDistance;    ///< Max allowed distance between split position and sliding linear fit position at the split x coordinate
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TrackSplittingTool::Factory::CreateAlgorithmTool() const
{
    return new TrackSplittingTool();
}

} // namespace lar_content

#endif // #ifndef TRACK_SPLITTING_TOOL_H
