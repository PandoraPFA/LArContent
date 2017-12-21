/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H
#define LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  CosmicRaySplittingAlgorithm class
 */
class CosmicRaySplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRaySplittingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

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
     *  @brief  Find the position of greatest scatter along a sliding linear fit
     * 
     *  @param  slidingFitResult the input sliding linear fit result
     *  @param  splitPosition the position of greatest scatter
     *  @param  splitDirection1 the direction vector just above the scatter position
     *  @param  splitDirection2 the direction vector just below the scatter position
     */
    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, pandora::CartesianVector &splitPosition,
        pandora::CartesianVector &splitDirection1, pandora::CartesianVector &splitDirection2) const;

    /**
     *  @brief  Find a second replacement cluster that aligns with the scatter of the first branch cluster
     * 
     *  @param  branchSlidingFitResult the sliding fit result for the branch cluster
     *  @param  replacementSlidingFitResult the sliding fit result for the replacement cluster
     *  @param  splitPosition the candidate split position on the branch cluster
     *  @param  splitDirection1 the first track direction just above the split position
     *  @param  splitDirection2 the second track direction just below the split position
     *  @param  lengthSquared1 figure of merit for association with first track direction
     *  @param  lengthSquared2 figure of merit for association with second track direction
     */
    pandora::StatusCode ConfirmSplitPosition(const TwoDSlidingFitResult &branchSlidingFitResult,
        const TwoDSlidingFitResult &replacementSlidingFitResult, const pandora::CartesianVector &splitPosition,
        const pandora::CartesianVector &splitDirection1, const pandora::CartesianVector &splitDirection2,
        float &lengthSquared1, float &lengthSquared2) const;

    /**
     *  @brief  Split a branch cluster for case of one replacement cluster
     * 
     *  @param  pBranchCluster the branch cluster
     *  @param  pReplacementCluster the replacement cluster
     *  @param  splitPosition the split position
     *  @param  forwardDirection the direction of the branch cluster just above the split position
     *  @param  backwardDirection the direction of the branch cluster just below the split position
     */
    pandora::StatusCode PerformSingleSplit(const pandora::Cluster *const pBranchCluster, const pandora::Cluster *const pReplacementCluster,
        const pandora::CartesianVector &splitPosition, const  pandora::CartesianVector &forwardDirection,
        const pandora::CartesianVector &backwardDirection) const;

    /**
     *  @brief  Split a branch cluster for case of two replacement clusters
     *
     *  @param  pBranchCluster the branch cluster
     *  @param  pReplacementCluster1 the first replacement cluster
     *  @param  pReplacementCluster2 the second replacement cluster
     *  @param  splitPosition the split position
     *  @param  splitDirection1 the direction of the branch cluster just above the split position
     *  @param  splitDirection2 the direction of the branch cluster just below the split position
     */
    pandora::StatusCode PerformDoubleSplit(const pandora::Cluster *const pBranchCluster, const pandora::Cluster *const pReplacementCluster1,
        const pandora::Cluster *const pReplacementCluster2, const pandora::CartesianVector &splitPosition,
        const pandora::CartesianVector &splitDirection1, const pandora::CartesianVector &splitDirection2) const;

    /**
     *  @brief  Get list of calo hits to move in order to split a branch cluster into two segments for case of one replacement cluster
     *
     *  @param  pBranchCluster the branch cluster
     *  @param  pReplacementCluster the replacement cluster
     *  @param  splitPosition the split position
     *  @param  forwardDirection the direction of the branch cluster just above the split position
     *  @param  backwardDirection the direction of the branch cluster just below the split position
     *  @param  caloHitList the output hits to be removed from the branch cluster
     */
    void GetCaloHitListToMove(const pandora::Cluster *const pBranchCluster, const pandora::Cluster *const pReplacementCluster,
        const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &forwardDirection,
        const pandora::CartesianVector &backwardDirection, pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Get lists of calo hits to move in order to split a branch cluster into two segments for case of two replacement clusters
     *
     *  @param  pBranchCluster the branch cluster
     *  @param  splitPosition the split position
     *  @param  splitDirection1 the direction of the branch cluster just above the split position
     *  @param  splitDirection2 the direction of the branch cluster just below the split position
     *  @param  caloHitList1 the first segment to be split from the branch cluster
     *  @param  caloHitList2 the second segment to be split from the branch cluster
     */
    void GetCaloHitListsToMove(const pandora::Cluster *const pBranchCluster, const pandora::CartesianVector &splitPosition,
        const pandora::CartesianVector &splitDirection1, const pandora::CartesianVector &splitDirection2,
        pandora::CaloHitList &caloHitList1, pandora::CaloHitList &caloHitList2) const;

    /**
     *  @brief  Identify crossed tracks formed from the branch cluster and its replacement cluster
     *
     *  @param  pBranchCluster the branch cluster
     *  @param  pReplacementCluster1 the first replacement cluster
     *  @param  pReplacementCluster2 the second replacement cluster
     *  @param  splitPosition the split position
     */
    bool IdentifyCrossedTracks(const pandora::Cluster *const pBranchCluster, const pandora::Cluster *const pReplacementCluster1,
        const pandora::Cluster *const pReplacementCluster2, const pandora::CartesianVector &splitPosition) const;

    /**
     *  @brief  Split the branch cluster and add hits to the replacement cluster
     *
     *  @param  pBranchCluster the branch cluster
     *  @param  caloHitListToMove the list of hits to be removed from the branch cluster and added to the replacement cluster
     *  @param  caloHitListToKeep to receive the list of calo hits to keep
     */
    pandora::StatusCode GetCaloHitListToKeep(const pandora::Cluster *const pBranchCluster, const pandora::CaloHitList &caloHitListToMove,
        pandora::CaloHitList &caloHitListToKeep) const;

    /**
     *  @brief  Split the branch cluster and add hits to the replacement cluster
     *
     *  @param  pBranchCluster the branch cluster
     *  @param  pReplacementCluster the replacement cluster
     *  @param  caloHitListToMove the list of hits to be removed from the branch cluster and added to the replacement cluster
     */
    pandora::StatusCode SplitCluster(const pandora::Cluster *const pBranchCluster, const pandora::Cluster *const pReplacementCluster,
        const pandora::CaloHitList &caloHitListToMove) const;

    float          m_clusterMinLength;                     ///< minimum length of clusters for this algorithm
    unsigned int   m_halfWindowLayers;                     ///< number of layers to use for half-window of sliding fit
    float          m_samplingPitch;                        ///< sampling pitch for walking along sliding linear fit
    float          m_maxCosSplittingAngle;                 ///< smallest scatter angle allowed when splitting cluster
    float          m_minCosMergingAngle;                   ///< largest relative angle allowed when merging clusters
    float          m_maxTransverseDisplacement;            ///< maximum transverse displacement of associated clusters
    float          m_maxLongitudinalDisplacement;          ///< maximum longitudinal displacement of associated clusters
    float          m_maxLongitudinalDisplacementSquared;   ///< longitudinal displacement squared
};

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H
