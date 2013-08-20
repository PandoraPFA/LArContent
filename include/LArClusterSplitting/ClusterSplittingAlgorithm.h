/**
 *  @file   LArContent/include/LArClusterSplitting/ClusterSplittingAlgorithm.h
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

    typedef std::list<pandora::Cluster*> ClusterSplittingList;

    /**
     *  @brief  Split cluster into its two fragments
     * 
     *  @param  pCluster address of the cluster
     *  @param  splitLayer the layer at which to perform the split
     *  @param  clusterSplittingList to receive the two cluster fragments
     */
    virtual pandora::StatusCode SplitCluster(pandora::Cluster *const pCluster, const unsigned int splitLayer, ClusterSplittingList &clusterSplittingList) const;

    /**
     *  @brief  Find position at which cluster should be split
     * 
     *  @param  pCluster address of the cluster
     *  @param  splitPosition the layer at which to perform the split
     */
    virtual pandora::StatusCode FindBestSplitPosition(const pandora::Cluster *const pCluster, pandora::CartesianVector &splitPosition) const = 0; 

   /**
     *  @brief  Whether a cluster is a split candidate
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    virtual bool IsPossibleSplit(const pandora::Cluster *const pCluster) const = 0;
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
