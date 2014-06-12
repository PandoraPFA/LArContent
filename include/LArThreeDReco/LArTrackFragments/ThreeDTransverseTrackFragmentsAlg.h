/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTrackFragmentsAlg.h
 *
 *  @brief  Header file for the three dimensional transverse track fragments algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRANSVERSE_TRACK_FRAGMENTS_ALGORITHM_H
#define LAR_THREE_D_TRANSVERSE_TRACK_FRAGMENTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArThreeDReco/LArTrackFragments/ThreeDFragmentsBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDTransverseTrackFragmentsAlg class
 */
class ThreeDTransverseTrackFragmentsAlg : public ThreeDFragmentsBaseAlgorithm
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

    float               m_minXOverlap;                      ///< requirement on minimum X overlap for associated clusters
    float               m_minXOverlapFraction;              ///< requirement on minimum X overlap fraction for associated clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTransverseTrackFragmentsAlg::Factory::CreateAlgorithm() const
{
    return new ThreeDTransverseTrackFragmentsAlg();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRANSVERSE_TRACK_FRAGMENTS_ALGORITHM_H
