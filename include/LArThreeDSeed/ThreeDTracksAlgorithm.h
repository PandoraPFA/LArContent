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
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);
    bool ExamineTensor();

    /**
     *  @brief  Calculate overlap result for special case with clusters at constant x
     * 
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     */
    void CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
        const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_constantXTreatment;       ///< Whether to use alternative OverlapResult calculation for constant x clusters
    float           m_pseudoChi2Cut;            ///< The pseudo chi2 cut to identify matched sampling points
    float           m_minMatchedFraction;       ///< The minimum matched sampling fraction to allow particle creation
    unsigned int    m_minMatchedPoints;         ///< The minimum number of matched sampling points to allow particle creation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACKS_ALGORITHM_H
