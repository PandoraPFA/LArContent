/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDLongitudinalTracksAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional longitudinal tracks algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_LONGITUDINAL_TRACKS_ALGORITHM_H
#define LAR_THREE_D_LONGITUDINAL_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDLongitudinalTracksAlgorithm class
 */
class ThreeDLongitudinalTracksAlgorithm : public ThreeDBaseAlgorithm<TrackOverlapResult>
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

    typedef LArClusterHelper::TwoDSlidingFitResult TwoDSlidingFitResult;
    typedef TwoDSlidingFitResult::LayerFitResultMap LayerFitResultMap;

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

    bool ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_vertexChi2Cut;                ///< 
    float           m_reducedChi2Cut;               ///< 
    float           m_cosOpeningAngleCut;           ///< 
    float           m_samplingPitch;                ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDLongitudinalTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDLongitudinalTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_LONGITUDINAL_TRACKS_ALGORITHM_H
