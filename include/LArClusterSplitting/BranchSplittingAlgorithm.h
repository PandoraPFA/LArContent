/**
 *  @file   LArContent/include/LArClusterSplitting/BranchSplittingAlgorithm.h
 * 
 *  @brief  Header file for the delta ray splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
#define LAR_BRANCH_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterSplitting/ClusterSplittingAndExtensionAlgorithm.h"

namespace lar
{

/**
 *  @brief  BranchSplittingAlgorithm class
 */
class BranchSplittingAlgorithm : public ClusterSplittingAndExtensionAlgorithm
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
    void FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, 
        const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector &branchSplitPosition, 
        pandora::CartesianVector &branchSplitDirection) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);   

    float         m_maxTransverseDisplacement;      ///< 
    float         m_maxLongitudinalDisplacement;    ///< 
    float         m_minCosRelativeAngle;            ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BranchSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BranchSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
