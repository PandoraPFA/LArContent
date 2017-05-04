/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/OvershootSplittingAlgorithm.h
 *
 *  @brief  Header file for the overshoot splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_OVERSHOOT_SPLITTING_ALGORITHM_H
#define LAR_OVERSHOOT_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitMultiSplitAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  OvershootSplittingAlgorithm class
 */
class OvershootSplittingAlgorithm : public TwoDSlidingFitMultiSplitAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    OvershootSplittingAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FindBestSplitPositions(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterSplittingMap) const;

    typedef std::pair<float, pandora::CartesianVector> MyTrajectoryPoint;
    typedef std::vector<MyTrajectoryPoint> MyTrajectoryPointList;

    /**
     *  @brief  Use sliding fit results to calculate intersections of clusters
     *
     *  @param  slidingFitResultMap the sliding fit result map
     *  @param  clusterIntersectionMap the map of cluster intersection points
     */
    void BuildIntersectionMap(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterIntersectionMap) const;

    /**
     *  @brief  Use intersection points to decide on splitting points
     *
     *  @param  slidingFitResultMap the sliding fit result map
     *  @param  clusterIntersectionMap the input map of cluster intersection points
     *  @param  sortedIntersectionMap the output map of sorted cluster intersection points
     */
    void BuildSortedIntersectionMap(const TwoDSlidingFitResultMap &slidingFitResultMap, const ClusterPositionMap &clusterIntersectionMap,
        ClusterPositionMap &sortedIntersectionMap) const;

    /**
     *  @brief  Select split positions from sorted list of candidate positions
     *
     *  @param  sortedIntersectionMap the input map of candidate split positions
     *  @param  clusterSplittingMap the output map of selected split positions
     */
    void PopulateSplitPositionMap(const ClusterPositionMap &sortedIntersectionMap, ClusterPositionMap &clusterSplittingMap) const;

    /**
     *  @brief  Sort pfos by number of constituent hits
     *
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByHitProjection(const MyTrajectoryPoint &lhs, const MyTrajectoryPoint &rhs);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_minClusterLength;          ///<
    float           m_maxClusterSeparation;      ///<
    float           m_minVertexDisplacement;     ///<
    float           m_maxIntersectDisplacement;
    float           m_minSplitDisplacement;
};

} // namespace lar_content

#endif // #ifndef LAR_OVERSHOOT_SPLITTING_ALGORITHM_H
