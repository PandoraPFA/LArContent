/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDTransverseTracksAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional transverse tracks algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H
#define LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

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

    /**
     *  @brief  ParticleComponent class
     */
    class ParticleComponent
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pClusterU the u cluster
         *  @param  pClusterV the v cluster
         *  @param  pClusterW the w cluster
         *  @param  trackOverlapResult the track overlap result
         */
        ParticleComponent(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW, const TrackOverlapResult &trackOverlapResult);

        /**
         *  @brief  Get the u cluster
         * 
         *  @return the u cluster
         */
        pandora::Cluster *GetClusterU() const;

        /**
         *  @brief  Get the v cluster
         * 
         *  @return the v cluster
         */
        pandora::Cluster *GetClusterV() const;

        /**
         *  @brief  Get the w cluster
         * 
         *  @return the w cluster
         */
        pandora::Cluster *GetClusterW() const;

        /**
         *  @brief  Get the track overlap result
         * 
         *  @return the track overlap result
         */
        const TrackOverlapResult &GetTrackOverlapResult() const;

    private:
        pandora::Cluster   *m_pClusterU;            ///< 
        pandora::Cluster   *m_pClusterV;            ///< 
        pandora::Cluster   *m_pClusterW;            ///< 
        TrackOverlapResult  m_trackOverlapResult;   ///< 
    };

    typedef std::vector<ParticleComponent> ParticleComponentList;

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

    typedef std::vector<TrackOverlapResult> TrackOverlapResultVector;

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

    /**
     *  @brief  Whether a requested (adjacent) element of the fit segment tensor exists
     * 
     *  @param  fitSegmentTensor the fit segment tensor
     *  @param  indexU the u index
     *  @param  indexV the v index
     *  @param  indexW the w index
     *  @param  incrementU whether to increment the u index
     *  @param  incrementV whether to increment the v index
     *  @param  incrementW whether to increment the w index
     *  @param  newIndexU to receive the new u index
     *  @param  newIndexV to receive the new v index
     *  @param  newIndexW to receive the new w index
     *  @param  trackOverlapResult to receive the track overlap result
     * 
     *  @return boolean
     */
    bool IsPresent(const FitSegmentTensor &fitSegmentTensor, const unsigned int indexU, const unsigned int indexV, const unsigned int indexW,
        const bool incrementU, const bool incrementV, const bool incrementW, unsigned int &newIndexU, unsigned int &newIndexV,
        unsigned int &newIndexW, TrackOverlapResult &trackOverlapResult) const;

    bool ExamineTensor();

    /**
     *  @brief  Build proto particle, starting with provided component and picking up any matched components in the overlap tensor
     * 
     *  @param  firstComponent the first particle component
     *  @param  protoParticle to receive the populated proto particle
     */
    void BuildProtoParticle(const ParticleComponent &firstComponent, ProtoParticle &protoParticle) const;

    /**
     *  @brief  Whether two particle components match, representing the same particle
     * 
     *  @param  firstComponent the first particle component
     *  @param  secondComponent the second particle component
     * 
     *  @return boolean
     */
    bool IsParticleMatch(const ParticleComponent &firstComponent, const ParticleComponent &secondComponent) const;

    /**
     *  @brief  Whether two clusters might match and represent the same particle
     * 
     *  @param  pFirstCluster the address of the first cluster
     *  @param  pSecondCluster the address of the second cluster
     * 
     *  @return boolean
     */
    bool IsPossibleMatch(pandora::Cluster *const pFirstCluster, pandora::Cluster *const pSecondCluster) const;

    /**
     *  @brief  Whether two clusters match, representing the same particle
     * 
     *  @param  pFirstCluster the address of the first cluster
     *  @param  pSecondCluster the address of the second cluster
     * 
     *  @return boolean
     */
    bool IsParticleMatch(pandora::Cluster *const pFirstCluster, pandora::Cluster *const pSecondCluster) const;

    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<pandora::Cluster*, LArClusterHelper::TwoDSlidingFitResult> SlidingFitResultMap;

    float               m_pseudoChi2Cut;            ///< The pseudo chi2 cut to identify matched sampling points
    float               m_minOverallMatchedFraction;///< The minimum matched sampling fraction to allow particle creation
    float               m_minSegmentMatchedFraction;///< The minimum segment matched sampling fraction to allow segment grouping
    unsigned int        m_minSegmentMatchedPoints;  ///< The minimum number of matched segment sampling points to allow segment grouping
    SlidingFitResultMap m_slidingFitResultMap;      ///< The sliding fit result map
};

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

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDTransverseTracksAlgorithm::ParticleComponent::ParticleComponent(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
        pandora::Cluster *pClusterW, const TrackOverlapResult &trackOverlapResult) :
    m_pClusterU(pClusterU),
    m_pClusterV(pClusterV),
    m_pClusterW(pClusterW),
    m_trackOverlapResult(trackOverlapResult)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ThreeDTransverseTracksAlgorithm::ParticleComponent::GetClusterU() const
{
    return m_pClusterU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ThreeDTransverseTracksAlgorithm::ParticleComponent::GetClusterV() const
{
    return m_pClusterV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ThreeDTransverseTracksAlgorithm::ParticleComponent::GetClusterW() const
{
    return m_pClusterW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackOverlapResult &ThreeDTransverseTracksAlgorithm::ParticleComponent::GetTrackOverlapResult() const
{
    return m_trackOverlapResult;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRANSVERSE_TRACKS_ALGORITHM_H
