/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDLongitudinalTracksAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional tracks algorithm class.
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
   

    typedef LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap LayerFitResultMap;

   

   
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

   

  
    void CalculateOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
        const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW,
	const pandora::CartesianVector vtxMerged3D, const pandora::CartesianVector endMerged3D, TrackOverlapResult& overlapResult);



    bool ExamineTensor();



    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDLongitudinalTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDLongitudinalTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACKS_ALGORITHM_H
