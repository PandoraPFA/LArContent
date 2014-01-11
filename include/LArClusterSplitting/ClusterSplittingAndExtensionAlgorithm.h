/**
 *  @file   LArContent/include/LArClusterSplitting/ClusterSplittingAndExtensionAlgorithm.h
 * 
 *  @brief  Header file for the branch splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_SPLITTING_AND_EXTENSION_ALGORITHM_H
#define LAR_CLUSTER_SPLITTING_AND_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  ClusterSplittingAndExtensionAlgorithm class
 */
class ClusterSplittingAndExtensionAlgorithm : public pandora::Algorithm
{

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

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
  
    /**
     *  @brief  Output the best split positions in branch and replacement clusters
     * 
     *  @param  branchSlidingFit the inputted sliding fit result for possible branch cluster
     *  @param  pReplacementCluster the inputted sliding fit result for possible replacement cluster
     *  @param  branchSplitPosition the outputted start position of the branch
     *  @param  branchSplitDirection the outputted start direction of the branch
     */
    virtual void FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, 
        const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector &branchSplitPosition, 
        pandora::CartesianVector &branchSplitDirection) const = 0;

private:
    
    /**
     *  @brief  Calculate RMS deviation of branch hits relative to the split direction
     * 
     *  @param  pCluster the input branch cluster
     *  @param  splitPosition the start position of the branch 
     *  @param  splitDirection the start direction of the branch
     */
    float CalculateBranchChi2(const pandora::Cluster* const pCluster, const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &splitDirection) const;

    /**
     *  @brief  Separate cluster into the branch hits to be split from the primary cluster
     * 
     *  @param  pCluster the input branch cluster
     *  @param  splitPosition the start position of the branch 
     *  @param  splitDirection the start direction of the branch
     *  @param  principalCaloHitList the hits to be added to the principal cluster
     *  @param  branchCaloHitList the hits to be split off into the output branch cluster
     */
    void SplitBranchCluster(const pandora::Cluster* const pCluster, const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &splitDirection, 
        pandora::CaloHitList &principalCaloHitList, pandora::CaloHitList &branchCaloHitList) const;

    unsigned int  m_shortHalfWindowLayers;          ///<
    unsigned int  m_longHalfWindowLayers;           ///< 
    float         m_minClusterLength;               ///< 
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_SPLITTING_AND_EXTENSION_ALGORITHM_H
