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

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  TwoDSlidingFitSplittingAndSwitchingAlgorithm class
 */
class TwoDSlidingFitSplittingAndSwitchingAlgorithm : public pandora::Algorithm
{

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);



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



    ////////////////


    /**
     *  @brief  TODO
     *
     *  @param  pCluster
     *  @param  splitPosition
     *  @param  splitDirection
     *  @param  firstCaloHitList the hits to be added to the first new cluster
     *  @param  secondCaloHitList the hits to be added to the second new cluster
     */
    void SplitCluster(pandora::Cluster *const pCluster,
        const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &splitDirection,
        pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;





    /**
     *  @brief  TODO
     *
     *  @param  pCluster1
     *  @param  pCluster2
     *  @param  splitPosition
     *  @param  firstDirection
     *  @param  secondDirection
     */
    pandora::StatusCode ReplaceClusters(pandora::Cluster *const pCluster1, pandora::Cluster *const pCluster2,
        const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &firstDirection,
        const pandora::CartesianVector &secondDirection) const;



    unsigned int  m_halfWindowLayers;               ///<


    float         m_minClusterLength;               ///<

};

} // namespace lar

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SWITCHING_ALGORITHM_H
