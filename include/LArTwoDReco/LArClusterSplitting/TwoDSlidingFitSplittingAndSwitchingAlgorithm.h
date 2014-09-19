/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting and switching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SWITCHING_ALGORITHM_H
#define LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SWITCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  TwoDSlidingFitSplittingAndSwitchingAlgorithm class
 */
class TwoDSlidingFitSplittingAndSwitchingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDSlidingFitSplittingAndSwitchingAlgorithm();

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Find the best split position and direction for a pair of clusters
     *
     *  @param slidingFit1 the sliding linear fit to the first cluster
     *  @param slidingFit2 the sliding linear fit to the second cluster
     *  @param splitPosition the output split position
     *  @param direction1 the output direction of the first new cluster
     *  @param direction2 the output direction of the second new cluster
     */
    virtual pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFit1, const TwoDSlidingFitResult &slidingFit2,
        pandora::CartesianVector &splitPosition, pandora::CartesianVector &direction1, pandora::CartesianVector &direction2) const = 0;

private:
    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Build the map of sliding fit results
     *
     *  @param  clusterVector the input cluster vector
     *  @param  slidingFitResultMap the output sliding fit result map
     */
    void BuildSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief  Split cluster at a given position and direction
     *
     *  @param  pCluster the cluster
     *  @param  splitPosition the position at which to split the cluster
     *  @param  splitDirection the direction of the un-crossed cluster
     *  @param  firstCaloHitList the hits to be added to the first new cluster
     *  @param  secondCaloHitList the hits to be added to the second new cluster
     */
    void SplitCluster(pandora::Cluster *const pCluster,
        const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &splitDirection,
        pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    /**
     *  @brief  Replace crossed clusters with un-crossed clusters
     *
     *  @param  pCluster1 the first cluster to be deleted
     *  @param  pCluster2 the second cluster to be deleted
     *  @param  splitPosition the split position
     *  @param  firstDirection the direction of the first new cluster
     *  @param  secondDirection the direction of the second new cluster
     */
    pandora::StatusCode ReplaceClusters(pandora::Cluster *const pCluster1, pandora::Cluster *const pCluster2,
        const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &firstDirection,
        const pandora::CartesianVector &secondDirection) const;

    unsigned int  m_halfWindowLayers;     ///< half window layers for sliding linear fot
    float         m_minClusterLength;     ///< minimum length of clusters
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SWITCHING_ALGORITHM_H
