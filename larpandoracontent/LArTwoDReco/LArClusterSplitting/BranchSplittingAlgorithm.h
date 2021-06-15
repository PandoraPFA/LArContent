/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h
 *
 *  @brief  Header file for the branch splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
#define LAR_BRANCH_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSplicingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  BranchSplittingAlgorithm class
 */
class BranchSplittingAlgorithm : public TwoDSlidingFitSplittingAndSplicingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    BranchSplittingAlgorithm();

private:
    void FindBestSplitPosition(const TwoDSlidingFitResult &branchSlidingFit, const TwoDSlidingFitResult &replacementSlidingFit,
        pandora::CartesianVector &replacementStartPosition, pandora::CartesianVector &branchSplitPosition,
        pandora::CartesianVector &branchSplitDirection) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxTransverseDisplacement;   ///<
    float m_maxLongitudinalDisplacement; ///<
    float m_minLongitudinalExtension;    ///<
    float m_minCosRelativeAngle;         ///<
    float m_projectionAngularAllowance;  ///<
};

} // namespace lar_content

#endif // #ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
