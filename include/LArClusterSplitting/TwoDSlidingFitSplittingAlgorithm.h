/**
 *  @file   LArContent/include/ClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h
 *
 *  @brief  Header file for the kink splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_ALGORITHM_H
#define LAR_TWO_D_SLIDING_FIT_SPLITTING_ALGORITHM_H 1

#include "LArClusterSplitting/ClusterSplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  TwoDSlidingFitSplittingAlgorithm class
 */
class TwoDSlidingFitSplittingAlgorithm : public ClusterSplittingAlgorithm
{
protected:
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Use sliding linear fit to identify the best split position
     *
     *  @param  slidingFitResult the input sliding fit result
     *  @param  splitPosition the best split position
     *
     *  @return pandora::StatusCode
     */
    virtual pandora::StatusCode FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult,
	pandora::CartesianVector& splitPosition) const = 0;

private:
    pandora::StatusCode SplitCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
	pandora::CaloHitList &secondCaloHitList) const;

    /**
     *  @brief  Use sliding linear fit to separate cluster into two fragments
     *
     *  @param  slidingFitResult the input sliding fit result
     *  @param  splitPosition the split position
     *  @param  firstCaloHitList the hits in the first cluster fragment
     *  @param  secondCaloHitList the hits in the second cluster fragment
     *
     *  @return pandora::StatusCode
     */
    pandora::StatusCode SplitCluster(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult,
	const pandora::CartesianVector& splitPosition, pandora::CaloHitList &firstCaloHitList,
	pandora::CaloHitList &secondCaloHitList) const;

    unsigned int    m_minClusterLayers;       ///<
    unsigned int    m_slidingFitHalfWindow;   ///<
};

} // namespace lar

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_ALGORITHM_H
