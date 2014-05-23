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
#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/ThreeDBaseAlgorithm.h"

namespace lar
{

class TensorManipulationTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDTransverseTracksAlgorithm class
 */
class ThreeDTransverseTracksAlgorithm : public ThreeDBaseAlgorithm<TransverseOverlapResult>
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
     *  @brief  Sort tensor elements by number of matched sampling points, using matched fraction then xoverlap span to resolve ties
     * 
     *  @param  lhs the first tensor element
     *  @param  rhs the second tensor element
     * 
     *  @return boolean
     */
    static bool SortByNMatchedSamplingPoints(const TensorType::Element &lhs, const TensorType::Element &rhs);

    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingFitResult &GetCachedSlidingFitResult(pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the layer window for the sliding linear fits
     * 
     *  @return the layer window for the sliding linear fits
     */
    unsigned int GetSlidingFitWindow() const;

    virtual void UpdateForNewCluster(pandora::Cluster *const pNewCluster);
    virtual void UpdateUponDeletion(pandora::Cluster *const pDeletedCluster);

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
    typedef std::map<unsigned int, TransverseOverlapResult> FitSegmentToOverlapResultMap;
    typedef std::map<unsigned int, FitSegmentToOverlapResultMap> FitSegmentMatrix;
    typedef std::map<unsigned int, FitSegmentMatrix> FitSegmentTensor;

    void PreparationStep();
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
    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int                m_nMaxTensorToolRepeats;    ///< The maximum number of repeat loops over tensor tools
    unsigned int                m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    float                       m_pseudoChi2Cut;            ///< The pseudo chi2 cut to identify matched sampling points
    float                       m_minSegmentMatchedFraction;///< The minimum segment matched sampling fraction to allow segment grouping
    unsigned int                m_minSegmentMatchedPoints;  ///< The minimum number of matched segment sampling points to allow segment grouping
    float                       m_minOverallMatchedFraction;///< The minimum matched sampling fraction to allow particle creation
    unsigned int                m_minOverallMatchedPoints;  ///< The minimum number of matched segment sampling points to allow particle creation
    float                       m_minSamplingPointsPerLayer;///< The minimum number of sampling points per layer to allow particle creation
    TwoDSlidingFitResultMap     m_slidingFitResultMap;      ///< The sliding fit result map

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

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ThreeDTransverseTracksAlgorithm::GetSlidingFitWindow() const
{
    return m_slidingFitWindow;
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
