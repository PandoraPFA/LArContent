/**
 *  @file   LArContent/include/LArClusterSplitting/BranchSplittingAlgorithm.h
 * 
 *  @brief  Header file for the branch splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
#define LAR_BRANCH_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  BranchSplittingAlgorithm class
 */
class BranchSplittingAlgorithm : public pandora::Algorithm
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

protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Populate cluster vector with good quality clusters
     * 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Remove a branch from a cluster and replace it with a second cluster
     * 
     *  @param  pBranchCluster the cluster containing a branch to be removed
     *  @param  pReplacementCluster the replacement cluster
     *  @param  branchStartPosition the position of the start of the branch
     *  @param  replacementStartPosition the position the start of the replacement
     */
    pandora::StatusCode ReplaceBranch(pandora::Cluster *const pBranchCluster, pandora::Cluster *const pReplacementCluster,
        const pandora::CartesianVector &branchStartPosition, const pandora::CartesianVector &replacementStartPosition) const;


private:
  
    /**
     *  @brief  Output the best split positions in branch and replacement clusters
     * 
     *  @param  branchSlidingFit the inputted sliding fit result for possible branch cluster
     *  @param  pReplacementCluster the inputted sliding fit result fot possible replacement cluster
     *  @param  branchStartPosition the outputted start position of the branch cluster
     *  @param  replacementStartPosition the outputted start position of the replacement cluster
     */
    pandora::StatusCode FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, 
        const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector& branchStartPosition, 
        pandora::CartesianVector& replacementStartPosition) const;


    unsigned int  m_halfWindowLayers;               ///< 
    unsigned int  m_stepSizeLayers;                 ///< 
    float         m_minClusterLength;               ///< 
    float         m_maxTransverseDisplacement;      ///< 
    float         m_maxLongitudinalDisplacement;    ///< 
    float         m_minCosRelativeAngle;            ///< 


};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BranchSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BranchSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
