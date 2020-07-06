/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewLongTracksTool.h
 *
 *  @brief  Header file for the long tracks tool class.
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_LONG_TRACKS_TOOL_H
#define TWO_VIEW_LONG_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewLongTracksTool class
 */
class TwoViewLongTracksTool : public TransverseMatrixTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewLongTracksTool();

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
    static bool IsLongerThanDirectConnections(IteratorList::const_iterator iIter, const MatrixType::ElementList &elementList,
        const unsigned int minMatchedSamplingPointRatio, const pandora::ClusterSet &usedClusters);

    bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find long tracks, hidden by simple ambiguities in the tensor
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindLongTracks(const MatrixType &overlapMatrix, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Select a list of long track-like elements from a set of connected tensor elements
     *
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectLongElements(const MatrixType::ElementList &elementList, const pandora::ClusterSet &usedClusters,
        IteratorList &iteratorList) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    float           m_minMatchingScore;                 ///< The min global matcing score for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_LONG_TRACKS_TOOL_H
