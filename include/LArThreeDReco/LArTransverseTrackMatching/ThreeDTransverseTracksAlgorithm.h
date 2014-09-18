/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h
 *
 *  @brief  Header file for the three dimensional transverse tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H
#define LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.h"

namespace lar_content
{

class TransverseTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDTransverseTracksAlgorithm class
 */
class ThreeDTransverseTracksAlgorithm : public ThreeDTracksBaseAlgorithm<TransverseOverlapResult>
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
     *  @brief  Default constructor
     */
    ThreeDTransverseTracksAlgorithm();

    /**
     *  @brief  Sort tensor elements by number of matched sampling points, using matched fraction then xoverlap span to resolve ties
     *
     *  @param  lhs the first tensor element
     *  @param  rhs the second tensor element
     *
     *  @return boolean
     */
    static bool SortByNMatchedSamplingPoints(const TensorType::Element &lhs, const TensorType::Element &rhs);

private:
    typedef std::map<unsigned int, TransverseOverlapResult> FitSegmentToOverlapResultMap;
    typedef std::map<unsigned int, FitSegmentToOverlapResultMap> FitSegmentMatrix;
    typedef std::map<unsigned int, FitSegmentMatrix> FitSegmentTensor;

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
        TransverseOverlapResult &overlapResult);

    /**
     *  @brief  Get the number of matched points for three fit segments and accompanying sliding fit results
     *
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     *  @param  fitSegmentTensor the fit segment tensor
     */
    void GetFitSegmentTensor(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
        const TwoDSlidingFitResult &slidingFitResultW, FitSegmentTensor &fitSegmentTensor) const;

    /**
     *  @brief  Get the overlap result for three fit segments and the accompanying sliding fit results
     *
     *  @param  fitSegmentU the fit segment u
     *  @param  fitSegmentV the fit segment v
     *  @param  fitSegmentW the fit segment w
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     *
     *  @return the overlap result
     */
    TransverseOverlapResult GetSegmentOverlap(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW,
        const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV, const TwoDSlidingFitResult &slidingFitResultW) const;

    /**
     *  @brief  Get the best overlap result, by examining the fit segment tensor
     *
     *  @param  fitSegmentTensor the fit segment tensor
     *
     *  @return the best overlap result
     */
    TransverseOverlapResult GetBestOverlapResult(const FitSegmentTensor &fitSegmentTensor) const;

    /**
     *  @brief  Get track overlap results for possible connected segments
     *
     *  @param  indexU the index u
     *  @param  indexV the index v
     *  @param  indexW the index w
     *  @param  transverseOverlapResultVector the transverse overlap result vector
     */
    void GetPreviousOverlapResults(const unsigned int indexU, const unsigned int indexV, const unsigned int indexW,
        FitSegmentTensor &fitSegmentSumTensor, TransverseOverlapResultVector &transverseOverlapResultVector) const;

    void ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<TransverseTensorTool*> TensorToolList;
    TensorToolList              m_algorithmToolList;        ///< The algorithm tool list

    unsigned int                m_nMaxTensorToolRepeats;    ///< The maximum number of repeat loops over tensor tools
    float                       m_pseudoChi2Cut;            ///< The pseudo chi2 cut to identify matched sampling points
    float                       m_minSegmentMatchedFraction;///< The minimum segment matched sampling fraction to allow segment grouping
    unsigned int                m_minSegmentMatchedPoints;  ///< The minimum number of matched segment sampling points to allow segment grouping
    float                       m_minOverallMatchedFraction;///< The minimum matched sampling fraction to allow particle creation
    unsigned int                m_minOverallMatchedPoints;  ///< The minimum number of matched segment sampling points to allow particle creation
    float                       m_minSamplingPointsPerLayer;///< The minimum number of sampling points per layer to allow particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseTensorTool class
 */
class TransverseTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDTransverseTracksAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTransverseTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTransverseTracksAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H
