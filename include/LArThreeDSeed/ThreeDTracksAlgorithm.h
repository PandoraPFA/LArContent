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

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDTracksAlgorithm class
 */
class ThreeDTracksAlgorithm : public ThreeDBaseAlgorithm
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
     *  @brief  
     */
    class OverlapResult
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  nMatchedSamplingPoints
         *  @param  nSamplingPoints
         */
        OverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints);

        /**
         *  @brief  Get the number of matched sampling points
         *
         *  @return the number of matched sampling points
         */
        unsigned int GetNMatchedSamplingPoints() const;

        /**
         *  @brief  Get the number of sampling points
         *
         *  @return the number of sampling points
         */
        unsigned int GetNSamplingPoints() const;

        /**
         *  @brief  Get the fraction of sampling points resulting in a match
         *
         *  @return the fraction of sampling points resulting in a match
         */
        float GetMatchedFraction() const;

    private:
        unsigned int    m_nMatchedSamplingPoints;       ///< The number of matched sampling points
        unsigned int    m_nSamplingPoints;              ///< The number of sampling points
        float           m_matchedFraction;              ///< The fraction of sampling points resulting in a match
    };

    void InitializeTensor();
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

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
    void UpdateTensor();
    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    OverlapTensor<OverlapResult>    m_overlapTensor;    ///< The overlap tensor
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTracksAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDTracksAlgorithm::OverlapResult::OverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints) :
    m_nMatchedSamplingPoints(nMatchedSamplingPoints),
    m_nSamplingPoints(nSamplingPoints),
    m_matchedFraction(0.f)
{
    if ((0 == m_nSamplingPoints) || (m_nMatchedSamplingPoints > m_nSamplingPoints))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_matchedFraction = static_cast<float>(m_nMatchedSamplingPoints) / static_cast<float>(m_nSamplingPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ThreeDTracksAlgorithm::OverlapResult::GetNMatchedSamplingPoints() const
{
    return m_nMatchedSamplingPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ThreeDTracksAlgorithm::OverlapResult::GetNSamplingPoints() const
{
    return m_nSamplingPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDTracksAlgorithm::OverlapResult::GetMatchedFraction() const
{
    return m_matchedFraction;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACKS_ALGORITHM_H
