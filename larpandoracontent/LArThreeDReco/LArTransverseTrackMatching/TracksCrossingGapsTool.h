/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.h
 *
 *  @brief  Header file for the long tracks tool class.
 *
 *  $Log: $
 */
#ifndef TRACKS_CROSSING_GAPS_TOOL_H
#define TRACKS_CROSSING_GAPS_TOOL_H 1

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TracksCrossingGapsTool class
 */
class TracksCrossingGapsTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TracksCrossingGapsTool();

    bool Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find tracks crossing gaps, with unambiguous connection but poor overlap due to gaps
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindTracks(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Select a list of track-like elements crossing a gap in one or more views from a set of connected tensor elements
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the set of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectElements(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
        const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Calculate the effective overlap fractions given a set of clusters, taking gaps into account
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the connected tensor element
     *  @param  xOverlapFractionU to receive the effective overlap fraction in the u view
     *  @param  xOverlapFractionV to receive the effective overlap fraction in the v view
     *  @param  xOverlapFractionW to receive the effective overlap fraction in the w view
     */
    void CalculateEffectiveOverlapFractions(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element,
        float &xOverlapFractionU, float &xOverlapFractionV, float &xOverlapFractionW) const;

    /**
     *  @brief  Calculate the effective overlap span given a set of clusters, taking gaps into account
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the connected tensor element
     *  @param  xMinEffU to receive the effective min u coordinate
     *  @param  xMaxEffU to receive the effective max u coordinate
     *  @param  xMinEffV to receive the effective min v coordinate
     *  @param  xMaxEffV to receive the effective max v coordinate
     *  @param  xMinEffW to receive the effective min w coordinate
     *  @param  xMaxEffW to receive the effective max w coordinate
     */
    void CalculateEffectiveOverlapSpan(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element,
        float &xMinEffU, float &xMaxEffU, float &xMinEffV, float &xMaxEffV, float &xMinEffW, float &xMaxEffW) const;

    /**
     *  @brief  Check whether there is any gap in the three U-V-W clusters combination
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the connected tensor element
     *  @param  xSample the x sampling position
     *  @param  gapInU to receive whether there is a gap in the u view
     *  @param  gapInV to receive whether there is a gap in the v view
     *  @param  gapInW to receive whether there is a gap in the w view
     *
     *  @return boolean
     */
    bool PassesGapChecks(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element, const float xSample,
        bool &gapInU, bool &gapInV, bool &gapInW) const;

    /**
     *  @brief  Check individually each cluster where a gap might be present
     *
     *  @param  xSample, the x coordinate we are checking
     *  @param  slidingFitResult1 the sliding fit result for the cluster in view 1
     *  @param  slidingFitResult2 the sliding fit result for the cluster in view 2
     *  @param  slidingFitResult3 the sliding fit result for the cluster in view 3
     *  @param  gapIn1 whether there is a gap in view 1
     *  @param  gapIn2 whether there is a gap in view 2
     *  @param  gapIn3 whether there is a gap in view 3
     *
     *  @return boolean
     */
    bool CheckXPositionInGap(const float xSample, const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
        const TwoDSlidingFitResult &slidingFitResult3, bool &gapIn1, bool &gapIn2, bool &gapIn3) const;

    /**
     *  @brief  Check whether a x position is at the end of the cluster
     *
     *  @param  xSample, the x coordinate of the point tested
     *  @param  pCluster the cluster we are interrogating for its extreme coordinates
     *
     *  @return boolean
     */
    bool IsEndOfCluster(const float xSample, const TwoDSlidingFitResult &slidingFitResult) const;

    float m_minMatchedFraction;                  ///< The min matched sampling point fraction for particle creation
    unsigned int m_minMatchedSamplingPoints;     ///< The min number of matched sampling points for particle creation
    float m_minXOverlapFraction;                 ///< The min x overlap fraction (in each view) for particle creation
    unsigned int m_minMatchedSamplingPointRatio; ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
    float m_maxGapTolerance;                     ///< The max gap tolerance
    float m_sampleStepSize;                      ///< The sampling step size used in association checks, units cm
    unsigned int m_maxAngleRatio;                ///< The max ratio allowed in the angle
};

} // namespace lar_content

#endif // #ifndef LONG_TRACKS_TOOL_H
