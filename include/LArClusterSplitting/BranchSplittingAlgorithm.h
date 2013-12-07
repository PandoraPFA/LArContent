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



    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;


    void FindBestBranchSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector& branchStartPosition, pandora::CartesianVector& replacementStartPosition) const;

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
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BranchSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BranchSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
