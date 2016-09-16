/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitMultiSplitAlgorithm.h
 *
 *  @brief  Header file for the 2D sliding fit multi-split algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_MULTI_SPLIT_ALGORITHM_H
#define LAR_TWO_D_SLIDING_FIT_MULTI_SPLIT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  TwoDTrackSplittingAlgorithm class
 */
class TwoDSlidingFitMultiSplitAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDSlidingFitMultiSplitAlgorithm();

protected:
    typedef std::unordered_map<const pandora::Cluster*, pandora::CartesianPointVector> ClusterPositionMap;

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    virtual void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const = 0;

    /**
     *  @brief  Determine best split positions based on sliding fit result
     *
     *  @param  slidingFitResultMap mapping from clusters to sliding fit results
     *  @param  clusterSplittingMap mapping from clusters to split positions
     */
    virtual void FindBestSplitPositions(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterSplittingMap) const = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Build the map of sliding fit results
     *
     *  @param  clusterVector the vector of selected clusters
     *  @param  halfWindowLayers the half-window to use for the sliding fits
     *  @param  slidingFitResultMap the sliding fit result map
     */
    void BuildSlidingFitResultMap(const pandora::ClusterVector &clusterVector, const unsigned int halfWindowLayers,
        TwoDSlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief  Split clusters
     *
     *  @param  slidingFitResultMap mapping from clusters to sliding fit results
     *  @param  clusterSplittingMap mapping from clusters to split positions
     */
    pandora::StatusCode SplitClusters(const TwoDSlidingFitResultMap &slidingFitResultMap, const ClusterPositionMap &clusterSplittingMap) const;

    /**
     *  @brief  Split cluster
     *
     *  @param  slidingFitResult input sliding fit result
     *  @param  splitPositionList vector of split positions
     */
    pandora::StatusCode SplitCluster(const TwoDSlidingFitResult &slidingFitResult, const pandora::CartesianPointVector &splitPositionList) const;

    unsigned int    m_slidingFitHalfWindow;      ///<
    std::string     m_inputClusterList;          ///<
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_MULTI_SPLIT_ALGORITHM_H
