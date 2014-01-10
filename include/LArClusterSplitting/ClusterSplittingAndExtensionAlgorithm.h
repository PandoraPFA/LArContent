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
     *  @param  branchStartPosition the position of the start of the branch
     *  @param  replacementStartPosition the position the start of the replacement
     */
    pandora::StatusCode ReplaceBranch(pandora::Cluster *const pBranchCluster, pandora::Cluster *const pReplacementCluster,
        const pandora::CartesianVector &branchStartPosition, const pandora::CartesianVector &branchStartDirection) const;
  
    /**
     *  @brief  Output the best split positions in branch and replacement clusters
     * 
     *  @param  branchSlidingFit the inputted sliding fit result for possible branch cluster
     *  @param  pReplacementCluster the inputted sliding fit result fot possible replacement cluster
     *  @param  branchStartPosition the outputted start position of the branch cluster
     *  @param  replacementStartPosition the outputted start position of the replacement cluster
     */
    virtual pandora::StatusCode FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, 
        const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector &branchStartPosition, 
        pandora::CartesianVector &branchStartDirection) const = 0;

private:

    unsigned int  m_shortHalfWindowLayers;          ///<
    unsigned int  m_longHalfWindowLayers;           ///< 
    float         m_minClusterLength;               ///< 
   
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_SPLITTING_AND_EXTENSION_ALGORITHM_H
