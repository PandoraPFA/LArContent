/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeViewLongitudinalTracksAlgorithm.h
 *
 *  @brief  Header file for the three view longitudinal tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_LONGITUDINAL_TRACKS_ALGORITHM_H
#define LAR_THREE_VIEW_LONGITUDINAL_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewTrackMatchingAlgorithm.h"

namespace lar_content
{

class LongitudinalTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewLongitudinalTracksAlgorithm class
 */
class ThreeViewLongitudinalTracksAlgorithm : public ThreeViewTrackMatchingAlgorithm<LongitudinalOverlapResult>
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeViewLongitudinalTracksAlgorithm();

private:
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     */
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        LongitudinalOverlapResult &overlapResult);

    /**
     *  @brief  Calculate the overlap result for given 3D vertex and end positions
     *
     *  @param  slidingFitResultU the sliding fit result u
     *  @param  slidingFitResultV the sliding fit result v
     *  @param  slidingFitResultW the sliding fit result w
     *  @param  vtxMerged3D the 3D vertex position
     *  @param  endMerged3D the 3D end position
     *  @param  overlapResult to receive the overlap result
     */
    void CalculateOverlapResult(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
        const TwoDSlidingFitResult &slidingFitResultW, const pandora::CartesianVector &vtxMerged3D, const pandora::CartesianVector &endMerged3D,
        TrackOverlapResult &overlapResult) const;

    void ExamineTensor();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<LongitudinalTensorTool*> TensorToolVector;
    TensorToolVector    m_algorithmToolVector;              ///< The algorithm tool vector

    unsigned int        m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools
    float               m_vertexChi2Cut;                    ///< The maximum allowed chi2 for associating end points from three views
    float               m_reducedChi2Cut;                   ///< The maximum allowed chi2 for associating hit positions from three views
    float               m_samplingPitch;                    ///< Pitch used to generate sampling points along tracks
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LongitudinalTensorTool class
 */
class LongitudinalTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewLongitudinalTracksAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewLongitudinalTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_LONGITUDINAL_TRACKS_ALGORITHM_H
