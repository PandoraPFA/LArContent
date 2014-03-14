/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the cluster splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <list>

namespace lar
{

/**
 *  @brief  ClusterSplittingAlgorithm class
 */
class ClusterSplittingAlgorithm : public pandora::Algorithm
{
protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Split cluster into two fragments
     *
     *  @param  pCluster address of the cluster
     *  @param  firstCaloHitList the hits in the first fragment
     *  @param  secondCaloHitList the hits in the second fragment
     */
    virtual pandora::StatusCode SplitCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
        pandora::CaloHitList &secondCaloHitList) const = 0;

private:
    typedef std::list<pandora::Cluster*> ClusterSplittingList;

    /**
     *  @brief  Split cluster into two fragments
     *
     *  @param  pCluster address of the cluster
     *  @param  clusterSplittingList to receive the two cluster fragments
     */
    pandora::StatusCode SplitCluster(pandora::Cluster *const pCluster, ClusterSplittingList &clusterSplittingList) const;
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
