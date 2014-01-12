/**
 *  @file   LArContent/include/LArClusterSplitting/DeltaRaySplittingAlgorithm.h
 * 
 *  @brief  Header file for the delta ray splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_SPLITTING_ALGORITHM_H
#define LAR_DELTA_RAY_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterSplitting/ClusterSplittingAndExtensionAlgorithm.h"

namespace lar
{

/**
 *  @brief  DeltaRaySplittingAlgorithm class
 */
class DeltaRaySplittingAlgorithm : public ClusterSplittingAndExtensionAlgorithm
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

    float         m_stepSize;                       ///< 
    float         m_maxTransverseDisplacement;      ///< 
    float         m_maxLongitudinalDisplacement;    ///< 
    float         m_minCosRelativeAngle;            ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRaySplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRaySplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_DELTA_RAY_SPLITTING_ALGORITHM_H
