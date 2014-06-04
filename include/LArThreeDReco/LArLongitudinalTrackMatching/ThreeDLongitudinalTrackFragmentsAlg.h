/**
 *  @file   LArContent/include/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTrackFragmentsAlg.h
 *
 *  @brief  Header file for the three dimensional longitudinal track fragments algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_LONGITUDINAL_TRACK_FRAGMENTS_ALGORITHM_H
#define LAR_THREE_D_LONGITUDINAL_TRACK_FRAGMENTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDFragmentsBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDLongitudinalTrackFragmentsAlg class
 */
class ThreeDLongitudinalTrackFragmentsAlg : public ThreeDFragmentsBaseAlgorithm
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
    void GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        pandora::CartesianPointList &projectedPositions) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDLongitudinalTrackFragmentsAlg::Factory::CreateAlgorithm() const
{
    return new ThreeDLongitudinalTrackFragmentsAlg();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_LONGITUDINAL_TRACK_FRAGMENTS_ALGORITHM_H
