/**
 *  @file   LArContent/include/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h
 *
 *  @brief  Header file for the three dimensional longitudinal tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_LONGITUDINAL_TRACKS_ALGORITHM_H
#define LAR_THREE_D_LONGITUDINAL_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.h"

namespace lar_content
{

class LongitudinalTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDLongitudinalTracksAlgorithm class
 */
class ThreeDLongitudinalTracksAlgorithm : public ThreeDTracksBaseAlgorithm<LongitudinalOverlapResult>
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Sort tensor elements by chi-squared
     * 
     *  @param  lhs the first tensor element
     *  @param  rhs the second tensor element
     * 
     *  @return boolean
     */
    static bool SortByChiSquared(const TensorType::Element &lhs, const TensorType::Element &rhs);

private:
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     */
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW,
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

    unsigned int    m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools
    float           m_vertexChi2Cut;                    ///< The maximum allowed chi2 for associating end points from three views
    float           m_reducedChi2Cut;                   ///< The maximum allowed chi2 for associating hit positions from three views
    float           m_samplingPitch;                    ///< Pitch used to generate sampling points along tracks

    typedef std::vector<LongitudinalTensorTool*> TensorToolList;
    TensorToolList  m_algorithmToolList;                ///< The algorithm tool list
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LongitudinalTensorTool class
 */
class LongitudinalTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDLongitudinalTracksAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDLongitudinalTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDLongitudinalTracksAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_LONGITUDINAL_TRACKS_ALGORITHM_H
