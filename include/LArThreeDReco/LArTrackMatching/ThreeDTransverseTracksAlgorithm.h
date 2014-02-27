/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional transverse tracks algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H
#define LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"

#include "LArThreeDReco/ThreeDBaseAlgorithm.h"

namespace lar
{

class TensorManipulationTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDTransverseTracksAlgorithm class
 */
class ThreeDTransverseTracksAlgorithm : public ThreeDBaseAlgorithm<TrackOverlapResult>
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

    typedef std::map<pandora::Cluster*, LArClusterHelper::TwoDSlidingFitResult> SlidingFitResultMap;

private:
    /**
     *  @brief  SlidingFitDirection enum
     */
    enum SlidingFitDirection
    {
        POSITIVE_IN_X,
        NEGATIVE_IN_X,
        UNCHANGED_IN_X,
        UNKNOWN
    };

    typedef LArClusterHelper::TwoDSlidingFitResult TwoDSlidingFitResult;
    typedef TwoDSlidingFitResult::LayerFitResultMap LayerFitResultMap;

    /**
     *  @brief  FitSegment class
     */
    class FitSegment
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  twoDSlidingFitResult the two dimensional sliding fit result
         *  @param  startLayerIter the layer fit result map iterator for the start layer
         *  @param  secondLayerIter the layer fit result map iterator for the end layer
         */
        FitSegment(const TwoDSlidingFitResult &twoDSlidingFitResult, LayerFitResultMap::const_iterator startLayerIter,
            LayerFitResultMap::const_iterator endLayerIter);

        /**
         *  @brief  Get start layer
         *
         *  @return the start layer
         */
        int GetStartLayer() const;

        /**
         *  @brief  Get end layer
         *
         *  @return the end layer
         */
        int GetEndLayer() const;

        /**
         *  @brief  Get the minimum x value
         *
         *  @return the minimum x value
         */
        float GetMinX() const;

        /**
         *  @brief  Get the maximum x value
         *
         *  @return the maximum x value
         */
        float GetMaxX() const;

        /**
         *  @brief  Get the start value of the u, v or w coordinate
         *
         *  @return the start value of the u, v or w coordinate
         */
        float GetStartValue() const;

        /**
         *  @brief  Get the end value of the u, v or w coordinate
         *
         *  @return the end value of the u, v or w coordinate
         */
        float GetEndValue() const;

        /**
         *  @brief  Whether the x coordinate increases between the start and end layers
         *
         *  @return boolean
         */
        bool IsIncreasingX() const;

    private:
        int         m_startLayer;               ///< The start layer
        int         m_endLayer;                 ///< The end layer
        float       m_minX;                     ///< The minimum x value
        float       m_maxX;                     ///< The maximum x value
        float       m_startValue;               ///< The start value of the u, v or w coordinate
        float       m_endValue;                 ///< The end value of the u, v or w coordinate
        bool        m_isIncreasingX;            ///< Whether the x coordinate increases between the start and end layers
    };

    typedef std::vector<FitSegment> FitSegmentList;
    typedef std::map<unsigned int, TrackOverlapResult> FitSegmentToOverlapResultMap;
    typedef std::map<unsigned int, FitSegmentToOverlapResultMap> FitSegmentMatrix;
    typedef std::map<unsigned int, FitSegmentMatrix> FitSegmentTensor;

    void PreparationStep();
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Get the fit segment list for a given sliding fit result
     * 
     *  @param  slidingFitResult the sliding fit result
     *  @param  fitSegmentList to receive the fit segment list
     */
    void GetFitSegmentList(const TwoDSlidingFitResult &slidingFitResult, FitSegmentList &fitSegmentList) const;

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
    TrackOverlapResult GetSegmentOverlap(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW,
        const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV, const TwoDSlidingFitResult &slidingFitResultW) const;

    /**
     *  @brief  Get the best overlap result, by examining the fit segment tensor
     * 
     *  @param  fitSegmentTensor the fit segment tensor
     * 
     *  @return the best overlap result
     */
    TrackOverlapResult GetBestOverlapResult(FitSegmentTensor &fitSegmentTensor) const;

    /**
     *  @brief  Get track overlap results for possible connected segments
     * 
     *  @param  indexU the index u
     *  @param  indexV the index v
     *  @param  indexW the index w
     *  @param  maxIndexU the max index u
     *  @param  maxIndexV the max index v
     *  @param  maxIndexW the max index w
     *  @param  trackOverlapResultVector the track overlap result vector
     */
    void GetPreviousOverlapResults(const unsigned int indexU, const unsigned int indexV, const unsigned int indexW, const unsigned int maxIndexU,
        const unsigned int maxIndexV, const unsigned int maxIndexW, FitSegmentTensor &fitSegmentSumTensor, TrackOverlapResultVector &trackOverlapResultVector) const;

    bool ExamineTensor();
    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float                       m_pseudoChi2Cut;            ///< The pseudo chi2 cut to identify matched sampling points
    float                       m_minOverallMatchedFraction;///< The minimum matched sampling fraction to allow particle creation
    float                       m_minSegmentMatchedFraction;///< The minimum segment matched sampling fraction to allow segment grouping
    unsigned int                m_minSegmentMatchedPoints;  ///< The minimum number of matched segment sampling points to allow segment grouping
    SlidingFitResultMap         m_slidingFitResultMap;      ///< The sliding fit result map

    typedef std::vector<TensorManipulationTool*> TensorManipulationToolList;
    TensorManipulationToolList  m_algorithmToolList;        ///< The algorithm tool list
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TensorManipulationTool class
 */
class TensorManipulationTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDTransverseTracksAlgorithm::TensorType TensorType;
    typedef ThreeDTransverseTracksAlgorithm::SlidingFitResultMap SlidingFitResultMap;

    /**
     *  @brief  Run the algorithm tool
     * 
     *  @param  slidingFitResultMap the sliding fit result map
     *  @param  overlapTensor the track overlap result tensor
     *  @param  protoParticleVector the proto particle vector
     */
    virtual pandora::StatusCode Run(const SlidingFitResultMap &slidingFitResultMap, TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTransverseTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTransverseTracksAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline int ThreeDTransverseTracksAlgorithm::FitSegment::GetStartLayer() const
{
    return m_startLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int ThreeDTransverseTracksAlgorithm::FitSegment::GetEndLayer() const
{
    return m_endLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTransverseTracksAlgorithm::FitSegment::GetMinX() const
{
    return m_minX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTransverseTracksAlgorithm::FitSegment::GetMaxX() const
{
    return m_maxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTransverseTracksAlgorithm::FitSegment::GetStartValue() const
{
    return m_startValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTransverseTracksAlgorithm::FitSegment::GetEndValue() const
{
    return m_endValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDTransverseTracksAlgorithm::FitSegment::IsIncreasingX() const
{
    return m_isIncreasingX;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H
