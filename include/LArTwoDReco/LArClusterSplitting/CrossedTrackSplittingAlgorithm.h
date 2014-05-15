/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CrossedTrackSplittingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H
#define LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h"

namespace lar
{

/**
 *  @brief  CrossedTrackSplittingAlgorithm class
 */
class CrossedTrackSplittingAlgorithm : public TwoDSlidingFitSplittingAndSwitchingAlgorithm
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

    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFit1, const TwoDSlidingFitResult &slidingFit2,
        pandora::CartesianVector &splitPosition, pandora::CartesianVector &direction1, pandora::CartesianVector &direction2) const;




    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CrossedTrackSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CrossedTrackSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H
