/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h
 *
 *  @brief  Header file for the delta ray splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_SPLITTING_ALGORITHM_H
#define LAR_DELTA_RAY_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSplicingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRaySplittingAlgorithm class
 */
class DeltaRaySplittingAlgorithm : public TwoDSlidingFitSplittingAndSplicingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRaySplittingAlgorithm();

private:
    void FindBestSplitPosition(const TwoDSlidingFitResult &branchSlidingFit, const TwoDSlidingFitResult &replacementSlidingFit,
        pandora::CartesianVector &replacementStartPosition, pandora::CartesianVector &branchSplitPosition,
        pandora::CartesianVector &branchSplitDirection) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_stepSize;                    ///<
    float m_maxTransverseDisplacement;   ///<
    float m_maxLongitudinalDisplacement; ///<
    float m_minCosRelativeAngle;         ///<
};

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_SPLITTING_ALGORITHM_H
