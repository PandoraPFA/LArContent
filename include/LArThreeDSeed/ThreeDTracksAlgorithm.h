/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDTracksAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional tracksalgorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRACKS_ALGORITHM_H
#define LAR_THREE_D_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDTracksAlgorithm class
 */
class ThreeDTracksAlgorithm : public ThreeDBaseAlgorithm<TrackOverlapResult>
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

    typedef LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap LayerFitResultMap;

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
        FitSegment(const LArClusterHelper::TwoDSlidingFitResult &twoDSlidingFitResult, LayerFitResultMap::const_iterator startLayerIter,
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

    /**
     *  @brief  SegmentComparison class
     */
    class SegmentComparison
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  fitSegmentU fit segment u
         *  @param  fitSegmentV fit segment v
         *  @param  fitSegmentW fit segment w
         */
        SegmentComparison(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW);

        /**
         *  @brief  Get fit segment u
         *
         *  @return the address of fit segment u
         */
        const FitSegment &GetFitSegmentU() const;

        /**
         *  @brief  Get fit segment v
         *
         *  @return the address of fit segment v
         */
        const FitSegment &GetFitSegmentV() const;

        /**
         *  @brief  Get fit segment w
         *
         *  @return the address of fit segment w
         */
        const FitSegment &GetFitSegmentW() const;

    private:
        FitSegment  m_fitSegmentU;              ///< Fit segment u
        FitSegment  m_fitSegmentV;              ///< Fit segment v
        FitSegment  m_fitSegmentW;              ///< Fit segment w
    };

    typedef std::vector<SegmentComparison> SegmentComparisonList;

    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Get the fit segment list for a given sliding fit result
     * 
     *  @param  slidingFitResult the sliding fit result
     *  @param  fitSegmentList to receive the fit segment list
     */
    void GetFitSegmentList(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, FitSegmentList &fitSegmentList) const;

    /**
     *  @brief  Get the segment comparison list for a give set of fit segment lists
     * 
     *  @param  fitSegmentListU the fit segment list u
     *  @param  fitSegmentListV the fit segment list v
     *  @param  fitSegmentListW the fit segment list w
     *  @param  segmentComparisonList to receive the segment comparison list
     */
    void GetSegmentComparisonList(const FitSegmentList &fitSegmentListU, const FitSegmentList &fitSegmentListV, const FitSegmentList &fitSegmentListW,
       SegmentComparisonList &segmentComparisonList) const;

    /**
     *  @brief  Get the overlap result for a given segment comparison objects and the accompanying sliding fit results
     * 
     *  @param  segmentComparison the segment comparison details
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     * 
     *  @return the overlap result
     */
    TrackOverlapResult CalculateOverlapResult(const SegmentComparison &segmentComparison, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
        const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW) const;

    /**
     *  @brief  Calculate overlap result for special case with clusters at constant x
     * 
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     */
    void CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
        const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW);

    bool ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_pseudoChi2Cut;            ///< The pseudo chi2 cut to identify matched sampling points
    float           m_minMatchedFraction;       ///< The minimum matched sampling fraction to allow particle creation
    unsigned int    m_minMatchedPoints;         ///< The minimum number of matched sampling points to allow particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTracksAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline int ThreeDTracksAlgorithm::FitSegment::GetStartLayer() const
{
    return m_startLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int ThreeDTracksAlgorithm::FitSegment::GetEndLayer() const
{
    return m_endLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTracksAlgorithm::FitSegment::GetMinX() const
{
    return m_minX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTracksAlgorithm::FitSegment::GetMaxX() const
{
    return m_maxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTracksAlgorithm::FitSegment::GetStartValue() const
{
    return m_startValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTracksAlgorithm::FitSegment::GetEndValue() const
{
    return m_endValue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDTracksAlgorithm::FitSegment::IsIncreasingX() const
{
    return m_isIncreasingX;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDTracksAlgorithm::SegmentComparison::SegmentComparison(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW) :
    m_fitSegmentU(fitSegmentU),
    m_fitSegmentV(fitSegmentV),
    m_fitSegmentW(fitSegmentW)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDTracksAlgorithm::FitSegment &ThreeDTracksAlgorithm::SegmentComparison::GetFitSegmentU() const
{
    return m_fitSegmentU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDTracksAlgorithm::FitSegment &ThreeDTracksAlgorithm::SegmentComparison::GetFitSegmentV() const
{
    return m_fitSegmentV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDTracksAlgorithm::FitSegment &ThreeDTracksAlgorithm::SegmentComparison::GetFitSegmentW() const
{
    return m_fitSegmentW;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACKS_ALGORITHM_H
