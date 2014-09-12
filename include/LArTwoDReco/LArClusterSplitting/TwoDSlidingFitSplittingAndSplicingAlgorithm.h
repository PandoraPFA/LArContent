/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSplicingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting and splicing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SPLICING_ALGORITHM_H
#define LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SPLICING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar_content
{

/**
 *  @brief  TwoDSlidingFitSplittingAndSplicingAlgorithm class
 */
class TwoDSlidingFitSplittingAndSplicingAlgorithm : public pandora::Algorithm
{

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  ClusterExtension class
     */
    class ClusterExtension
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pBranchCluster  the address of the branch cluster
         *  @param  pReplacementCluster  the address of the replacement cluster
         *  @param  replacementVertex  the start position of the replacement cluster
         *  @param  branchVertex  the split position of the branch cluster
         *  @param  branchDirection  the split direction of the branch cluster
         */
        ClusterExtension(pandora::Cluster *const pBranchCluster, pandora::Cluster *const pReplacementCluster,
            const pandora::CartesianVector &replacementVertex, const pandora::CartesianVector &branchVertex,
            const pandora::CartesianVector &branchDirection);

        /**
         *  @brief  return the address of the branch Cluster
         */
        pandora::Cluster* GetBranchCluster() const;

        /**
         *  @brief  return the address of the replacement Cluster
         */
        pandora::Cluster* GetReplacementCluster() const;

        /**
         *  @brief  return the start position of the replacement cluster
         */
        const pandora::CartesianVector &GetReplacementVertex() const;

        /**
         *  @brief  return the split position of the branch cluster
         */
        const pandora::CartesianVector &GetBranchVertex() const;

        /**
         *  @brief  return the split direction of the branch cluster
         */
        const pandora::CartesianVector &GetBranchDirection() const;

        private:
        pandora::Cluster*           m_pBranchCluster;            ///<
        pandora::Cluster*           m_pReplacementCluster;       ///<
        pandora::CartesianVector    m_replacementVertex;         ///<
        pandora::CartesianVector    m_branchVertex;              ///<
        pandora::CartesianVector    m_branchDirection;           ///<
    };

    typedef std::vector<ClusterExtension> ClusterExtensionList;

    /**
     *  @brief  Output the best split positions in branch and replacement clusters
     *
     *  @param  branchSlidingFit the inputted sliding fit result for possible branch cluster
     *  @param  pReplacementCluster the inputted sliding fit result for possible replacement cluster
     *  @param  replacementStartPosition the outputted start position of the replacement
     *  @param  branchSplitPosition the outputted start position of the branch
     *  @param  branchSplitDirection the outputted start direction of the branch
     */
    virtual void FindBestSplitPosition(const TwoDSlidingFitResult &branchSlidingFit, 
        const TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector &replacementStartPosition,
        pandora::CartesianVector &branchSplitPosition, pandora::CartesianVector &branchSplitDirection) const = 0;

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
     *  @param  halfWindowLayers the half-window to use for the sliding fits
     *  @param  slidingFitResultMap the output sliding fit result map
     */
    void BuildSlidingFitResultMap(const pandora::ClusterVector &clusterVector, const unsigned int halfWindowLayers,
        TwoDSlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief  Build a list of candidate splits
     *
     *  @param  clusterVector the input cluster vector
     *  @param  branchResultMap the sliding fit result map for branch clusters
     *  @param  replacementResultMap the sliding fit result map for replacement clusters
     *  @param  clusterExtensionList the output list of candidate splits
     */
    void BuildClusterExtensionList(const pandora::ClusterVector &clusterVector, const TwoDSlidingFitResultMap &branchResultMap,
        const TwoDSlidingFitResultMap &replacementResultMap, ClusterExtensionList &clusterExtensionList) const;

    /**
     *  @brief  Finalize the list of candidate splits
     *
     *  @param  inputList the input list of possible splits
     *  @param  branchResultMap the sliding fit result map for branch clusters
     *  @param  replacementResultMap the sliding fit result map for replacement clusters
     *  @param  outputList the output list of definite splits
     */
    void PruneClusterExtensionList(const ClusterExtensionList &inputList, const TwoDSlidingFitResultMap &branchResultMap,
        const TwoDSlidingFitResultMap &replacementResultMap, ClusterExtensionList &outputList) const;

    /**
     *  @brief  Calculate RMS deviation of branch hits relative to the split direction
     *
     *  @param  pCluster the input branch cluster
     *  @param  splitPosition the start position of the branch
     *  @param  splitDirection the start direction of the branch
     */
    float CalculateBranchChi2(const pandora::Cluster* const pCluster, const pandora::CartesianVector &splitPosition, 
        const pandora::CartesianVector &splitDirection) const;

    /**
     *  @brief  Separate cluster into the branch hits to be split from the primary cluster
     *
     *  @param  pCluster the input branch cluster
     *  @param  splitPosition the start position of the branch
     *  @param  splitDirection the start direction of the branch
     *  @param  principalCaloHitList the hits to be added to the principal cluster
     *  @param  branchCaloHitList the hits to be split off into the output branch cluster
     */
    void SplitBranchCluster(const pandora::Cluster* const pCluster, const pandora::CartesianVector &splitPosition, 
        const pandora::CartesianVector &splitDirection, pandora::CaloHitList &principalCaloHitList, 
        pandora::CaloHitList &branchCaloHitList) const;

    /**
     *  @brief  Run the machinary that performs the cluster splitting and extending
     *
     *  @param  splitList the input list of candidate splits
     *  @param  branchResultMap the sliding fit result map for branch clusters
     *  @param  replacementResultMap the sliding fit result map for replacement clusters
     */
    pandora::StatusCode RunSplitAndExtension(const ClusterExtensionList &splitList, TwoDSlidingFitResultMap &branchResultMap,
        TwoDSlidingFitResultMap &replacementResultMap) const;

   /**
     *  @brief  Remove a branch from a cluster and replace it with a second cluster
     *
     *  @param  pBranchCluster the cluster containing a branch to be removed
     *  @param  pReplacementCluster the replacement cluster
     *  @param  branchSplitPosition the position at the start of the branch
     *  @param  branchSplitDirection the direction at the start of the branch
     */
    pandora::StatusCode ReplaceBranch(pandora::Cluster *const pBranchCluster, pandora::Cluster *const pReplacementCluster,
        const pandora::CartesianVector &branchSplitPosition, const pandora::CartesianVector &branchSplitDirection) const;

    unsigned int  m_shortHalfWindowLayers;          ///<
    unsigned int  m_longHalfWindowLayers;           ///<
    float         m_vetoDisplacement;               ///<
    float         m_minClusterLength;               ///<
    bool          m_runCosmicMode;                  ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoDSlidingFitSplittingAndSplicingAlgorithm::ClusterExtension::ClusterExtension(pandora::Cluster *const pBranchCluster,
    pandora::Cluster *const pReplacementCluster, const pandora::CartesianVector &replacementVertex,
    const pandora::CartesianVector &branchVertex, const pandora::CartesianVector &branchDirection) :
    m_pBranchCluster(pBranchCluster),
    m_pReplacementCluster(pReplacementCluster),
    m_replacementVertex(replacementVertex),
    m_branchVertex(branchVertex),
    m_branchDirection(branchDirection)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster* TwoDSlidingFitSplittingAndSplicingAlgorithm::ClusterExtension::GetBranchCluster() const
{
    return m_pBranchCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster* TwoDSlidingFitSplittingAndSplicingAlgorithm::ClusterExtension::GetReplacementCluster() const
{
    return m_pReplacementCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TwoDSlidingFitSplittingAndSplicingAlgorithm::ClusterExtension::GetReplacementVertex() const
{
    return m_replacementVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TwoDSlidingFitSplittingAndSplicingAlgorithm::ClusterExtension::GetBranchVertex() const
{
    return m_branchVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TwoDSlidingFitSplittingAndSplicingAlgorithm::ClusterExtension::GetBranchDirection() const
{
    return m_branchDirection;
}

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_SPLITTING_AND_SPLICING_ALGORITHM_H
