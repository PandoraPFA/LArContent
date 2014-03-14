/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/SeedBranchGrowingAlgorithm.h
 * 
 *  @brief  Header file for the seed branch growing algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_BRANCH_GROWING_ALGORITHM_H
#define LAR_SEED_BRANCH_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "SeedGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  SeedBranchGrowingAlgorithm class
 */
class SeedBranchGrowingAlgorithm : public SeedGrowingAlgorithm
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
    virtual void GetCandidateClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    virtual AssociationType AreClustersAssociated(const pandora::Cluster *const pClusterSeed, const pandora::Cluster *const pCluster) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedBranchGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedBranchGrowingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_SEED_BRANCH_GROWING_ALGORITHM_H
